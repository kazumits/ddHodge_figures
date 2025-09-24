import sys, time, warnings, logging, json
import numpy as np
import anndata as ad

# suppress numba warnings
warnings.simplefilter('ignore')
import dynamo as dyn
from dynamo.vectorfield.utils import compute_curl, compute_divergence
from dynamo.simulation.ODE import two_genes_motif_jacobian

M = int(sys.argv[1])
N = int(sys.argv[2])

def f_jac(X):
    if X.ndim == 1: X = X.reshape((1, -1))
    J = np.zeros((X.shape[1], X.shape[1], X.shape[0]))
    for ind, i in enumerate(X):
        J[:, :, ind] = two_genes_motif_jacobian(i[0], i[1])
    return J

def gen_data(N=1000):
    #adata = dyn.sim.Simulator(motif="twogenes", cell_num=N)
    adata = ad.read_h5ad(f"datasets/twogenes_N{N}.h5ad")
    adata.obsm['X_pca'], adata.obsm['velocity_pca'] = adata.X, adata.layers['velocity']
    adata.var['use_for_dynamics'] = True # i.e., use all
    adata.uns['PCs'] = np.identity(2) # bypass PCA i.e., original space
    return adata

def gen_runner(M=100,silence=True):
    def runner(adata,silence=silence): 
        if(silence): logging.disable(logging.CRITICAL)
        dyn.vf.VectorField(
            adata, basis='pca', velocity_key='velocity',
            pot_curl_div=False, M=M, lstsq_method='scipy',
        )
        dyn.vf.jacobian(adata,basis='pca',store_in_adata=True)
        dyn.vf.curl(adata,basis='pca')
        dyn.vf.divergence(adata,basis='pca')
        if(silence): logging.disable(logging.NOTSET)
        return adata
    return runner

def eval_mses(adata):
    mse = lambda x, y: np.mean((x - y)**2)
    J_a = f_jac(adata.X)
    curl_a = compute_curl(f_jac, adata.X)
    div_a = compute_divergence(f_jac, adata.X, vectorize_size=1)
    detJ_a = [np.linalg.det(J_a[:,:,i]) for i in range(adata.n_obs)]
    return {
        'mse_div': mse(adata.obs.divergence_pca, div_a),
        'mse_curl': mse(adata.obs.curl_pca, curl_a),
        'mse_detJ': mse(adata.obs.jacobian_det_pca, detJ_a)
    }

runner = gen_runner(M)
adata = gen_data(N)
start = time.perf_counter()
adata = runner(adata)
end = time.perf_counter()

print(json.dumps({'M': M, 'N': N, 'time': end-start} | eval_mses(adata)))

