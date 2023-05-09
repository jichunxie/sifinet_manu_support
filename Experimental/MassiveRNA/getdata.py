import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp_sparse
import scipy.io as sio

### https://stackoverflow.com/questions/46733052/read-hdf5-file-into-numpy-array
### https://cf.10xgenomics.com/supp/cell-exp/megacell_tutorial.html


hf = h5py.File('1M_neurons_filtered_gene_bc_matrices_h5.h5', 'r')
indptr = hf.get('mm10/indptr')
data = hf.get('mm10/data')
indices = hf.get('mm10/indices')
shape = hf.get('mm10/shape')

matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
sio.mmwrite("1M_matrix.mtx", matrix)
