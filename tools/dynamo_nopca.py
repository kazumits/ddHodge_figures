import sys 
import anndata as ad
import dynamo as dyn
import numpy as np

adata = ad.read_h5ad(sys.argv[1])
adata.obsm['X_pca'], adata.obsm['velocity_pca'] = adata.X, adata.layers['velocity']
adata.var['use_for_dynamics'] = True # i.e., use all
adata.uns['PCs'] = np.identity(adata.X.shape[1]) # bypass PCA i.e., original space

dyn.vf.VectorField(
    adata, basis='pca', velocity_key='velocity', pot_curl_div=False,
    M=1000, lstsq_method='scipy'
)

dyn.vf.jacobian(adata, basis='pca', store_in_adata=True)
dyn.vf.curl(adata,basis='pca')
dyn.vf.divergence(adata,basis='pca')
adata.write(sys.argv[1])
print(adata)
