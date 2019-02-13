
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.manifold import TSNE

from scvi.dataset import AnnDataset
from scvi.dataset import CortexDataset, RetinaDataset, PbmcDataset
from scvi.models import *
from scvi.inference import UnsupervisedTrainer


# In[2]:


# create a python script from the notebook
get_ipython().system('jupyter nbconvert scVI_testing.ipynb --to script')


# In[3]:


#gene_dataset = AnnDataset(filename="subsampled_CLUESImmVar_nonorm.h5ad", save_path="./")
gene_dataset = PbmcDataset(save_path="../data")


# In[ ]:


n_epochs_all = None
#save_path = 'data/'
show_plot = True

n_epochs=400 if n_epochs_all is None else n_epochs_all
#n_epochs=4 if n_epochs_all is None else n_epochs_all
lr=1e-3
use_batches=True
use_cuda=True


# In[ ]:


print(gene_dataset.n_batches)


# In[ ]:


vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * use_batches)
trainer = UnsupervisedTrainer(vae,
                              gene_dataset,
                              train_size=0.75,
                              use_cuda=use_cuda,
                              frequency=5)
trainer.train(n_epochs=n_epochs, lr=lr)


# In[ ]:


n_samples_tsne = 1000
trainer.train_set.show_t_sne(n_samples=n_samples_tsne, color_by='labels')


# In[ ]:


# Plotting the likelihood change across the 50 epochs of training: blue for training error and orange for testing error. 

ll_train = trainer.history["ll_train_set"]
ll_test = trainer.history["ll_test_set"]
x = np.linspace(0,50,(len(ll_train)))
plt.plot(x, ll_train)
plt.plot(x, ll_test)
plt.ylim(min(ll_train)-50, 3500)


# In[ ]:


print("Entropy batch mixing :", trainer.train_set.entropy_batch_mixing())


# In[ ]:


# obtaining latent space in the same order as the input data
n_samples_tsne = 1000
trainer.train_set.show_t_sne(n_samples=n_samples_tsne, color_by='batches and labels', save_name='latent_proj.png')


# In[ ]:


latent_data, batch_indices, labels = trainer.train_set.get_latent()


# In[ ]:


pd.DataFrame(latent_data).to_csv("latent_proj.csv", sep=",", header=True, index=False)
covar_labels = pd.DataFrame(batch_indices, columns=["batch_cov"])
covar_labels["ct_cov"] = labels
covar_labels.to_csv("latent_proj_batch_labels.csv", sep=",", header=True, index=False)


# In[ ]:


get_ipython().run_cell_magic('R', '-i latent_data -i batch_indices -o batch_estimate', 'library(kBET)\n\nbatch_estimate <- kBET(latent_data, batch_indices, do.pca=FALSE, n_repeat=20, verbose=TRUE)')


# In[ ]:


batch_estimate

