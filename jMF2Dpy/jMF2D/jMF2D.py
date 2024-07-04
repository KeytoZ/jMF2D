import numpy as np
from .utils import *
import time

def prepare_for_Deconvolution(sc_data, sc_annotation, sp_data, sp_location=None,
                              filter=True, species="Human", min_cells=3, max_genes=5000,
                              select_hvg=True, top_n_genes=5000,
                              preprocess_=True, scale=True,
                              construct=True, n_neighbors=6):
    """
    pandas.DataFrame format
        sc_data: gene * cell
        sc_annotaion: cell * celltype, and aligned in the cells and sc_data
        sp_data: gene * spot
        sp_location: Contains the corresponding x and y coordinates for each spot

    parameters
        filter: Whether or not sc data is pre-filtered
        species: Human or Mouse to filter mitochondrial genes, min_cells and max_genes refer to scanpy
        select_hvg: Whether or not the HVG is calculated, top_n_genes: The number of HVG genes
        preprocess_: scanpy standard data processing process, normal_total, log, and scale
        conctruct: Whether or not construct the spatial network
        n_neighbors: The number of neighbors built in the spatial network

    return:
        X: The DataFrame of the spatial transcriptome with dimension g*spot
        Y: The DataFrame of the signature expression of celltype with dimension g*celltype
        L: The Laplacian matrix of topological structure of spatial position
    """
    if(filter):
        print("## filter anndata...")
        sc_adata = sc.AnnData(X=sc_data.T)
        print(sc_adata)
        sc_adata.obs['annotation'] = sc_annotation
        sc_adata = filter_sc(sc_adata, species=species, min_cells=min_cells, max_genes=max_genes)
        sc_data = sc_adata.to_df().T
        sc_annotation = sc_adata.obs['annotation'].values

    intersected_gene = list(set(sp_data.index.values) & set(sc_data.index.values))
    sc_data = sc_data.loc[intersected_gene, :]

    if(select_hvg):
        print("## select the high varibale genes...")
        sc_adata = sc.AnnData(X=sc_data.T)
        sc.pp.highly_variable_genes(sc_adata, flavor="seurat_v3", n_top_genes=top_n_genes)
        sc_adata = sc_adata[:, sc_adata.var['highly_variable']]
        sc_data = sc_adata.to_df().T

    if(preprocess_):
        print("## preprocess the scRNA-seq anndata...")
        sc_adata = sc.AnnData(X=sc_data.T)
        preprocess(sc_adata, scale)
        sc_data = sc_adata.to_df().T

    sc_data = sc_data.T
    sc_data['annotation'] = sc_annotation
    print("## Each celltype contains the number of cells: \n", sc_data.groupby('annotation').count())
    # The mean expression is calculated by celltype
    Y = sc_data.groupby('annotation').mean().T


    # sp_data
    sp_data = sp_data.loc[Y.index.values, :]

    if(preprocess_):
        print("## preprocess the spatial anndata...")
        sp_adata = sc.AnnData(X=sp_data.T)
        preprocess(sp_adata)
        sp_data = sp_adata.to_df().T
    # gene alignment
    X = sp_data.T
    X = X[Y.index.values]
    X = X.T

    if construct:
        print("## construct the spatial network...")
        # spatial loacation
        W = conctruct_KNN(sp_location, n_neighbors=n_neighbors)
        D = np.diag(W.sum(axis=1))
        L = D - W
    else:
        L = np.eye(X.shape[1])

    return X, Y, L

def runDeconvolution(X, Y, L2, alpha=1, beta=1, gamma=1, loss=1e-3, max_iter=500,
    seed=0, initial_B=1, initial_F=0, sym=0, normbyspot=0):
    """
    There are many parameters that can be fine-tuned, but the core is alpha,beta,gamma,
    and the meaning of each parameter is as follows
        Parameters:
            X : The data matrix of the spatial transcriptome with dimension g*spot
            Y: The signature expression of celltype with dimension g*celltype
            L2: The Laplacian matrix of topological structure of spatial position
            alpha: The L21 norm weight parameter
            beta: The constraint strength of cell similarity on base matrix B
            gamma: The constraint strength of space topology on coefficient matrix F
            loss: Iteration termination loss
            max_iter: The upper limit for the number of iterations
            seed: If set to 0, the random seed is not set, otherwise it is set
            initial_B: The basis matrix B is initialized in a different way, 0 or 1 or any other data represents a different way
            initial_F: The coefficient matrix F is initialized in a different way, 0 or any other number
            sym_Z: Ensure that Z is a symmetric matrix mode, and 0 or any other number represents a different mode
            norm_F: Whether to normalize F by spot during iteration, 0 is false and 1 is true
    """
    X = X.values
    Y = Y.values
    s = X.shape[1]
    g = X.shape[0]
    k = Y.shape[1]
    if seed!=0:
        np.random.seed(seed)
    if initial_B==0:
        B = np.random.rand(g,k)
    elif initial_B==1:
        B = Y
    else:
        B = Y + np.random.rand(g,k)
    if initial_F==0:
        F = np.random.rand(k,s)
    else:
        F = np.linalg.lstsq(B,X)[0]
    E = np.random.rand(k,s)
    Z = np.ones((k,k))
    C = np.ones((k,k))
    T1 = np.zeros((k,k))
    T2 = np.zeros((k,s))
    delta1 = 1
    delta2 = 1

    ## Update parameters iteratively
    converge = 0
    iter = 0
    start_time = time.time()
    while converge==0 and iter<max_iter:
        # Update Z
        Zk = Z.copy()
        Z = Z - np.diag(np.diag(Z))
        Z_num = Y.T @ Y + delta1 * C - T1
        Z_den = Y.T @ Y @ Z + delta1 * Z
        Z_den[Z_den == 0] = np.finfo(float).eps
        Z = Z * (Z_num / Z_den)
        Z = Z / Z.sum(axis=1)[:,None]
        Z = Z - np.diag(np.diag(Z)-1)
        if sym==0:
            Z = (Z + Z.conj().T) / 2
        else:
            Z = np.triu(Z, k=0)
            Z = Z + Z.conj().T - np.eye(k)
        tz = np.linalg.norm(Z-Zk,ord="fro")

        # cal L1
        L1 = np.diag(np.sum(Z,axis=0)) - Z

        # Update C
        Ck = C.copy()
        Crow = np.sqrt(np.sum(np.abs(np.conj(C.conj().T)) * C.conj().T, axis=0))
        D1 = np.diag((1/2)*Crow)
        C_num = delta1 * Z + T1
        C_den = alpha * D1 @ C + delta1 * C
        C_den[C_den == 0] = np.finfo(float).eps
        C = C * (C_num / C_den)
        tc = np.linalg.norm(C-Ck,ord="fro")

        # Update T1
        T1 = T1 + delta1 * (Z-C)

        # Update B
        Bk = B.copy()
        B_num = X @ F.T
        B_den = B @ F @ F.T + (beta / 2) * (B @ L1 + B @ L1.T)
        B_den[B_den == 0] = np.finfo(float).eps
        B = B * (B_num / B_den)
        tb = np.linalg.norm(B-Bk,ord="fro")
        B = np.maximum(B, 0)
        B[np.isnan(B)] = 0

        # Update F
        Fk = F.copy()
        F_num = B.T @ X + delta2 * E
        F_den = B.T @ B @ F + delta2 * F + T2
        F_den[F_den == 0] = np.finfo(float).eps
        F = F * (F_num / F_den)
        if normbyspot==1:
            F = F / F.sum(axis=0)[:, None]
        F = np.real(F)
        tf = np.linalg.norm(F-Fk,ord="fro")
        F = np.maximum(F, 0)

        # Update E
        Ek = E.copy()
        Erow = np.sqrt(np.sum(np.abs(E.conj().T) * E.T, axis=0))
        D2 = np.diag(1 / 2 * Erow)
        E_num = delta2 * F + T2
        E_den = (gamma / 2) * (E @ L2 + E @ L2.T) + alpha * D2 @ E + delta2 * E
        E_den[E_den == 0] = np.finfo(float).eps
        E = E * (E_num / E_den)
        E = np.maximum(E, 0)
        te = np.linalg.norm(E-Ek,ord="fro")

        # Update T2
        T2 = T2 + delta2 * (F-E)

        # cal loss
        tol_loss = max(np.linalg.norm(B-Bk,ord='fro'),tf) / max(np.linalg.norm(B,ord='fro'), np.linalg.norm(F,ord='fro'))
        if tol_loss < loss:
            converge=1
        else:
            iter = iter + 1
        print("iter:{} loss:{} ".format(iter,tol_loss))
    end_time = time.time()
    print("Run times (s):{:.2f}".format(end_time-start_time))
    return F, Z

