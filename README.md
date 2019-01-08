# Index_model
Computational model to calculate tumor index (TI) using gene expression profile of samples

Step 1: install R packages (GSA,randomForest, preprocessCore, pheatmap)
Step 2: run model_training.R to train model using data in input.csv and as input 
Only expression profiles of WT livers and tumours will be used to train model, tumour index WT liver was assumed as -1, tumour as 1
Step 3: run model_test.R to calculate tumour index of samples using their gene expression profiles
