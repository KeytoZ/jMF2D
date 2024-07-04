import pandas as pd
import numpy as np
import scipy.stats as st
import scanpy as sc
from scipy import stats
from sklearn.neighbors import NearestNeighbors
import ot

## Evaluate the performance of algorithm
def cal_ssim(im1, im2, M):
    """
        calculate the SSIM value between two arrays.

    Parameters
        -------
        im1: array1, shape dimension = 2
        im2: array2, shape dimension = 2
        M: the max value in [im1, im2]
    """

    assert len(im1.shape) == 2 and len(im2.shape) == 2
    assert im1.shape == im2.shape
    mu1 = im1.mean()
    mu2 = im2.mean()
    sigma1 = np.sqrt(((im1 - mu1) ** 2).mean())
    sigma2 = np.sqrt(((im2 - mu2) ** 2).mean())
    sigma12 = ((im1 - mu1) * (im2 - mu2)).mean()
    k1, k2, L = 0.01, 0.03, M
    C1 = (k1 * L) ** 2
    C2 = (k2 * L) ** 2
    C3 = C2 / 2
    l12 = (2 * mu1 * mu2 + C1) / (mu1 ** 2 + mu2 ** 2 + C1)
    c12 = (2 * sigma1 * sigma2 + C2) / (sigma1 ** 2 + sigma2 ** 2 + C2)
    s12 = (sigma12 + C3) / (sigma1 * sigma2 + C3)
    ssim = l12 * c12 * s12

    return ssim

def scale_max(df):
    """
        Divided by maximum value to scale the data between [0,1].
        Please note that these datafrmae are scaled data by column.

        Parameters
        -------
        df: dataframe, each col is a feature.

    """

    result = pd.DataFrame()
    for label, content in df.items():
        content = content / content.max()
        result = pd.concat([result, content], axis=1)
    return result

def scale_z_score(df):
    """
        scale the data by Z-score to conform the data to the standard normal distribution, that is, the mean value is 0, the standard deviation is 1, and the conversion function is 0.
        Please note that these datafrmae are scaled data by column.

        Parameters
        -------
        df: dataframe, each col is a feature.
    """

    result = pd.DataFrame()
    for label, content in df.items():
        content = stats.zscore(content)
        content = pd.DataFrame(content, columns=[label])
        result = pd.concat([result, content], axis=1)
    return result

def scale_plus(df):
    """
        Divided by the sum of the data to scale the data between (0,1), and the sum of data is 1.
        Please note that these datafrmae are scaled data by column.

        Parameters
        -------
        df: dataframe, each col is a feature.
    """

    result = pd.DataFrame()
    for label, content in df.items():
        content = content / content.sum()
        result = pd.concat([result, content], axis=1)
    return result

def ssim(raw, impute, scale='scale_max'):
    ## This was used for calculating the SSIM value between two arrays.

    if scale == 'scale_max':
        raw = scale_max(raw)
        impute = scale_max(impute)
    else:
        print('Please note you do not scale data by max')
    if raw.shape[1] == impute.shape[1]:
        result = pd.DataFrame()
        for label in raw.columns:
            raw_col = raw.loc[:, label]
            impute_col = impute.loc[:, label]

            M = [raw_col.max(), impute_col.max()][raw_col.max() > impute_col.max()]
            raw_col_2 = np.array(raw_col)
            raw_col_2 = raw_col_2.reshape(raw_col_2.shape[0], 1)

            impute_col_2 = np.array(impute_col)
            impute_col_2 = impute_col_2.reshape(impute_col_2.shape[0], 1)

            ssim = cal_ssim(raw_col_2, impute_col_2, M)

            ssim_df = pd.DataFrame(ssim, index=["SSIM"], columns=[label])
            result = pd.concat([result, ssim_df], axis=1)
        return result
    else:
        print("columns error")

def pearsonr(raw, impute, scale=None):
    ## This was used for calculating the PCC between two arrays.
    if raw.shape[1] == impute.shape[1]:
        result = pd.DataFrame()
        for label in raw.columns:
            raw_col = raw.loc[:, label]
            impute_col = impute.loc[:, label]
            pearsonr, _ = st.pearsonr(raw_col, impute_col)
            pearson_df = pd.DataFrame(pearsonr, index=["Pearson"], columns=[label])
            result = pd.concat([result, pearson_df], axis=1)
        return result

def JS(raw, impute, scale='scale_plus'):
    ## This was used for calculating the JS value between two arrays.

    if scale == 'scale_plus':
        raw = scale_plus(raw)
        impute = scale_plus(impute)
    else:
        print('Please note you do not scale data by plus')
    if raw.shape[1] == impute.shape[1]:
        result = pd.DataFrame()
        for label in raw.columns:
            raw_col = raw.loc[:, label]
            impute_col = impute.loc[:, label]

            M = (raw_col + impute_col) / 2
            KL = 0.5 * st.entropy(raw_col, M) + 0.5 * st.entropy(impute_col, M)
            KL_df = pd.DataFrame(KL, index=["JS"], columns=[label])

            result = pd.concat([result, KL_df], axis=1)
        return result

def RMSE(raw, impute, scale='zscore'): # 'zscore'
    ## This was used for calculating the RMSE value between two arrays.

    if scale == 'zscore':
        raw = scale_z_score(raw)
        impute = scale_z_score(impute)
    else:
        print('Please note you do not scale data by zscore')
    if raw.shape[1] == impute.shape[1]:
        result = pd.DataFrame()
        for label in raw.columns:
            raw_col = raw.loc[:, label]
            impute_col = impute.loc[:, label]
            RMSE = np.sqrt(((raw_col - impute_col) ** 2).mean())
            RMSE_df = pd.DataFrame(RMSE, index=["RMSE"], columns=[label])

            result = pd.concat([result, RMSE_df], axis=1)
        return result


## Prepare the datas for jMF2D to Deconvolution
def preprocess(adata, scale=True):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    if scale:
        sc.pp.scale(adata, zero_center=False, max_value=10)

def filter_sc(adata, species="Human", min_cells=3, max_genes=5000):
    if(max_genes > 0):
        sc.pp.filter_cells(adata, max_genes=max_genes)
    if(min_cells > 0):
        sc.pp.filter_genes(adata, min_cells=min_cells)
    if species=="Human":
        mit_gene = "MT-"
    else:
        mit_gene = "mt-"
    adata.var['mt'] = adata.var_names.str.startswith(mit_gene)
    adata = adata[:, ~adata.var['mt']]
    return adata

def conctruct_KNN(position, n_neighbors):
    if isinstance(position, pd.DataFrame):
        position = position.values
    n_spot = position.shape[0]
    print("spot num ", n_spot)

    nbrs = NearestNeighbors(algorithm='ball_tree').fit(position)
    graph_out = nbrs.kneighbors_graph(n_neighbors=n_neighbors,
                                      mode="distance")
    graph_out.data = 1 / graph_out.data
    # row_normalize
    for start_ptr, end_ptr in zip(graph_out.indptr[:-1], graph_out.indptr[1:]):
        row_sum = graph_out.data[start_ptr:end_ptr].sum()
        if row_sum != 0:
            graph_out.data[start_ptr:end_ptr] /= row_sum
        # print(f"normalized sum from ptr {start_ptr} to {end_ptr} "
        #       f"({end_ptr - start_ptr} entries)",
        #       np.sum(graph_out.data[start_ptr:end_ptr]))
    return graph_out.A

def refine_label(adata, radius=50, key='label'):
    n_neigh = radius
    new_type = []
    old_type = adata.obs[key].values

    # calculate distance
    position = adata.obsm['spatial']
    distance = ot.dist(position, position, metric='euclidean')
    n_cell = distance.shape[0]

    for i in range(n_cell):
        vec = distance[i, :]
        index = vec.argsort()
        neigh_type = []
        for j in range(1, n_neigh + 1):
            neigh_type.append(old_type[index[j]])
        max_type = max(neigh_type, key=neigh_type.count)
        new_type.append(max_type)

    new_type = [str(i) for i in list(new_type)]
    # adata.obs['label_refined'] = np.array(new_type)

    return new_type

def Moran_I(genes_exp,x, y, k=5, knn=True):
    XYmap=pd.DataFrame({"x": x, "y":y})
    if knn:
        XYnbrs = NearestNeighbors(n_neighbors=k, algorithm='auto',metric = 'euclidean').fit(XYmap)
        XYdistances, XYindices = XYnbrs.kneighbors(XYmap)
        W = np.zeros((genes_exp.shape[0],genes_exp.shape[0]))
        for i in range(0,genes_exp.shape[0]):
            W[i,XYindices[i,:]]=1
        for i in range(0,genes_exp.shape[0]):
            W[i,i]=0

    I = pd.Series(index=genes_exp.columns, dtype="float64")
    for k in genes_exp.columns:
        X_minus_mean = np.array(genes_exp[k] - np.mean(genes_exp[k]))
        X_minus_mean = np.reshape(X_minus_mean,(len(X_minus_mean),1))
        Nom = np.sum(np.multiply(W,np.matmul(X_minus_mean,X_minus_mean.T)))
        Den = np.sum(np.multiply(X_minus_mean,X_minus_mean))
        I[k] = (len(genes_exp[k])/np.sum(W))*(Nom/Den)
    return I

def Geary_C(genes_exp,x, y, k=5, knn=True):
    XYmap=pd.DataFrame({"x": x, "y":y})
    if knn:
        XYnbrs = NearestNeighbors(n_neighbors=k, algorithm='auto',metric = 'euclidean').fit(XYmap)
        XYdistances, XYindices = XYnbrs.kneighbors(XYmap)
        W = np.zeros((genes_exp.shape[0],genes_exp.shape[0]))
        for i in range(0,genes_exp.shape[0]):
            W[i,XYindices[i,:]]=1
        for i in range(0,genes_exp.shape[0]):
            W[i,i]=0

    C = pd.Series(index=genes_exp.columns, dtype="float64")
    for k in genes_exp.columns:
        X=np.array(genes_exp[k])
        X_minus_mean = X - np.mean(X)
        X_minus_mean = np.reshape(X_minus_mean,(len(X_minus_mean),1))
        Xij=np.array([X,]*X.shape[0]).transpose()-np.array([X,]*X.shape[0])
        Nom = np.sum(np.multiply(W,np.multiply(Xij,Xij)))
        Den = np.sum(np.multiply(X_minus_mean,X_minus_mean))
        C[k] = (len(genes_exp[k])/(2*np.sum(W)))*(Nom/Den)
    return C