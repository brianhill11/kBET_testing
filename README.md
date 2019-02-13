# kBET_testing
Code for testing kBET metric

## Creating the Anaconda environment 
To replicate the Anaconda environment, make sure conda is installed and then run this command:
```
conda env create -f environment.yml
```

## Generating data
The generate_data.ipynb notebook contains code to generate single-cell data using datasets included in the [scVI package](https://github.com/YosefLab/scVI). It creates data in 3 formats: 
1. cells x genes matrix (`"prefix".csv` file)
2. genes x cells matrix (`"prefix"_T.csv` file)
3. cells x genes h5ad file (`"prefix".h5ad`)

Additionally, a cells x covariates file `"prefix"_batch_labels.csv` is created, where the first column is the batch label for each cell, and the second column is the predicted cell-type label for each cell. 

Note that the code currently dumps _dense_ matrices, which means for large datasets the files can be big and it can take a while to dump the data to disk. Also note that these matrices compress fairly well, so a future optimization could be to write them as compressed files, assuming the methods can read compresssed files. 

## Running [Seurat](https://github.com/satijalab/seurat)
The R Seurat is installed using CRAN. 

## Running [DCA](https://github.com/theislab/dca) 
The DCA package is installed with conda. 

## Running [scVI](https://github.com/YosefLab/scVI)
The scVI package is installed with conda. 

## Running [kBET](https://github.com/theislab/kBET)
The kBET R package is installed from Github using the following commands: 
```
library(devtools)
install_github('theislab/kBET')
```
