from utils import *
import matlab.engine
import matlab
import scipy

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
    eng = matlab.engine.start_matlab()
    F, Z = eng.deconvolution(matlab.double(X), matlab.double(Y), matlab.double(L2), matlab.double([alpha]), matlab.double([beta]), matlab.double([gamma]), loss, max_iter,
    seed, initial_B, initial_F, sym, normbyspot, nargout=2)
    eng.quit()
    return F, Z

