#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd

import scvi
from scvi.dataset.pbmc import PbmcDataset
from scvi.dataset.loom import RetinaDataset

import anndata


# In[2]:


# create a python script from the notebook
get_ipython().system('jupyter nbconvert generate_data.ipynb --to script')


# ## Configurations

# In[3]:


# this sets the file name prefix for generated files
prefix = "retina_eq_ct"

# choose a dataset to use
#gene_dataset = PbmcDataset()
gene_dataset = RetinaDataset()

# optionally, balance the cell-type composition between the batches
equal_cell_types = False


# ## Get cell-type and batch labels from scVI dataset object
# For each dataset we need: 
#     (1) an gene expression count matrix 
#     (2) cell type labels
#     (3) batch labels

# In[4]:


# get cell type labels
cell_types = [gene_dataset.cell_types[i] for i in gene_dataset.labels.ravel()]

# labels = pd.DataFrame(gene_dataset.batch_indices, columns=["batch_cov"])
# labels["ct_cov"] = cell_types
# labels.to_csv(prefix + "_batch_labels.csv", sep=",", header=True, index=False)


# In[5]:


# batch labels
gene_dataset.batch_indices


# In[6]:


print(len(gene_dataset.batch_indices))


# ## Get the expression data, batch labels, and cell type labels into an AnnData object

# In[14]:


ad = anndata.AnnData(gene_dataset.X, obs={"batch_cov": list([i[0] for i in gene_dataset.batch_indices]), "ct_cov": cell_types})


# In[17]:


print(ad.X.shape)


# In[18]:


# for each batch, look at distribution of cell types
for batch in np.unique(ad.obs["batch_cov"]):
    print("="*30)
    print(batch)
    print("="*30)
    print(ad.obs[ad.obs["batch_cov"]==batch]["ct_cov"].value_counts())
    print(ad.obs[ad.obs["batch_cov"]==batch]["ct_cov"].value_counts()/ad.obs[ad.obs["batch_cov"]==batch].shape[0])


# ## Optional: subsample cells so that cell type compositions match between batches
# 

# In[19]:


if equal_cell_types:
    # subsample the cells so that cell type proportions match in each batch
    min_num_cells_per_cell_type = {}
    for ct in ad.obs["ct_cov"].unique():
        min_num_cells = np.min(ad.obs[ad.obs["ct_cov"]==ct]["batch_cov"].value_counts())
        min_num_cells_per_cell_type[ct] = min_num_cells
    print(min_num_cells_per_cell_type)

    print(ad.shape)

    # for a batch, if we find it has more than the minimum number of a particular cell type, downsample
    final_indices = []
    for batch in ad.obs["batch_cov"].unique():
        for ct in ad.obs["ct_cov"].unique():
    #        if sc_data[(sc_data.obs["batch_cov"]==batch) & (sc_data.obs["ct_cov"]==ct)].shape[0] > min_num_cells_per_cell_type[ct]:
            indices = ad[(ad.obs["batch_cov"]==batch) & (ad.obs["ct_cov"]==ct)].to_df().index
            indices_to_keep = np.random.choice(indices, size=min_num_cells_per_cell_type[ct], replace=False)
            final_indices = final_indices + list(indices_to_keep)

    print(len(final_indices))
    ad = ad[final_indices, :]
    print(ad.shape)
    
    # for each batch, look at distribution of cell types
    for batch in ad.obs["batch_cov"].unique():
        print("="*30)
        print(batch)
        print("="*30)
        print(ad.obs[ad.obs["batch_cov"]==batch]["ct_cov"].value_counts())
        print(ad.obs[ad.obs["batch_cov"]==batch]["ct_cov"].value_counts()/ad.obs[ad.obs["batch_cov"]==batch].shape[0])


# In[10]:


# create the cells x genes matrix csv files
try: 
    pd.DataFrame(ad.X.todense(), columns=ad.var_names).to_csv(prefix + ".csv", sep=",", header=True, index=False)
except AttributeError:
    pd.DataFrame(ad.X, columns=ad.var_names).to_csv(prefix + ".csv", sep=",", header=True, index=False)


# In[11]:


# this step is for generating data for DCA (requires genes x cells matrix)
try: 
    pd.DataFrame(ad.X.todense(), columns=ad.var_names).T.to_csv(prefix + "_T.csv", sep=",", header=True, index=False)
except AttributeError:
    pd.DataFrame(ad.X, columns=ad.var_names).T.to_csv(prefix + "_T.csv", sep=",", header=True, index=False)


# In[12]:


# get cell type labels from object
cell_types = ad.obs["ct_cov"]

#labels = pd.DataFrame(gene_dataset.batch_indices, columns=["batch_cov"])
labels = pd.DataFrame(ad.obs["batch_cov"], columns=["batch_cov"])
labels["ct_cov"] = cell_types
print(labels.head(20))
print(labels.shape)
labels.to_csv(prefix + "_batch_labels.csv", sep=",", header=True, index=False)


# In[13]:


# write AnnDataset to h5ad file
ad.write(prefix + ".h5ad")


# In[ ]:




