# Index_model
Computational model to calculate tumor index (TI) using gene expression profile of samples

Step 1: install R packages (GSA,randomForest, preprocessCore, pheatmap, glmnet)
Step 2: run model_training.R to train model using data in data.csv and as input 
Only expression profiles of WT livers and tumours will be used to train model, tumour index WT liver was assumed as -1, tumour as 1
Step 3: run model_testing.R to calculate tumour index of samples using their gene expression profiles
calculate index of all samples in data.csv
