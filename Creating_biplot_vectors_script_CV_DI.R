#Download the package from Bioconductor
#if (!requireNamespace('BiocManager', quietly = TRUE))
  #install.packages('BiocManager')

#BiocManager::install('PCAtools')
#Note: to install development version direct from GitHub:
  
 # if (!requireNamespace('devtools', quietly = TRUE))
    #install.packages('devtools')

#devtools::install_github('kevinblighe/PCAtools')
#Load the package into R session
library(PCAtools)


#set working directory where matrix and other files are located
  
  setwd("/Users/michelle/Documents/German Lab Research Mac/Cv_DI_Liver_Ls_PC/Cv_DI_Transcriptome/updated PCA plot Feb 2022")

    #load primary data from the matrix file!!
  
    primary_data = read.table("Cv_DI_FINAL_GENES_2022_redo.isoform.counts.matrix", header=T, com='', row.names=1, check.names=F, sep='\t')
  primary_data = as.matrix(primary_data)


  ##replace NA with 0 in data matrix
  primary_data[is.na(primary_data)] <- 0
  
  
  PCA_CV_DI <- pca(primary_data, removeVar = 0.1)
  
  #Change species names in the biplot commands and set Diet as a factor! 
  PCA_CV_DI[["yvars"]]
  Diet <- PCA_CV_DI[["yvars"]]
  Diet <- as.factor(Diet)
  
  screeplot(PCA_CV_DI, axisLabSize = 18, titleLabSize = 22)
  biplot(PCA_CV_DI, lab = Diet, gridlines.major = FALSE, gridlines.minor = FALSE, labSize = 2, pointSize = 10, shape = Diet, shapekey = c('Cv118_DI'= 19,'Cv119_DI'= 19, 'Cv305_DI'= 19, 'Cv317_DI'= 19, 'Cv306_DI'= 19,'Cv312_DI'= 19, 'Cv307_DI'= 19,'Cv313_DI_redo_RSEM'= 19), colkey = c('Cv118_DI'= 'gray50','Cv119_DI'= 'gray50', 'Cv305_DI'= 'plum3', 'Cv317_DI'= 'plum3', 'Cv306_DI'= 'steelblue','Cv312_DI'= 'steelblue', 'Cv307_DI'= 'tomato','Cv313_DI_redo_RSEM'= 'tomato'))
  #showloadings for vectors
  biplot(PCA_CV_DI, lab = Diet, showLoadings = TRUE,sizeLoadingsNames = 5, gridlines.major = FALSE, gridlines.minor = FALSE, labSize = 2, pointSize = 10, shape = Diet, shapekey = c('Cv118_DI'= 19,'Cv119_DI'= 19, 'Cv305_DI'= 19, 'Cv317_DI'= 19, 'Cv306_DI'= 19,'Cv312_DI'= 19, 'Cv307_DI'= 19,'Cv313_DI_redo_RSEM'= 19), colkey = c('Cv118_DI'= 'gray50','Cv119_DI'= 'gray50', 'Cv305_DI'= 'plum3', 'Cv317_DI'= 'plum3', 'Cv306_DI'= 'steelblue','Cv312_DI'= 'steelblue', 'Cv307_DI'= 'tomato','Cv313_DI_redo_RSEM'= 'tomato'))
  
  #clean up graph so labels and loading names are not there
  plot.margin=unit(c(10,10,10,10),"cm") 
  biplot(PCA_CV_DI, lab = Diet, labSize = 0, drawConnectors = FALSE, showLoadings = TRUE, showLoadingsNames = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, axisLabSize=30, borderWidth= 3, pointSize = 15, shape = Diet, shapekey = c('Cv118_DI'= 19,'Cv119_DI'= 19, 'Cv305_DI'= 19, 'Cv317_DI'= 19, 'Cv306_DI'= 19,'Cv312_DI'= 19, 'Cv307_DI'= 19,'Cv313_DI_redo_RSEM'= 19), colkey = c('Cv118_DI'= 'gray50','Cv119_DI'= 'gray50', 'Cv305_DI'= 'plum3', 'Cv317_DI'= 'plum3', 'Cv306_DI'= 'steelblue','Cv312_DI'= 'steelblue', 'Cv307_DI'= 'tomato','Cv313_DI_redo_RSEM'= 'tomato'))
  #could never figure out how to center figure so just manually labeled axis in word doc
  
 