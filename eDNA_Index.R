### Reference: Kelly, R. P., Shelton, A. O. & Gallego, R. Understanding PCR processes to draw meaningful conclusions from environmental DNA studies. Sci. Rep. 9, 12133–14 (2019).
### URL: https://github.com/invertdna/eDNA_Process_Simulations/blob/master/FunctionsForSimulationPaper.R

#function to calculate eDNA index
eDNAINDEX <- function(x) { #where x is a dataframe with taxa/OTUs/etc in rows, and samples in columns
  rowMax <- function(x){apply(x, MARGIN = 1, FUN = max)}
  temp <- sweep(x, MARGIN = 2, STATS = colSums(x), FUN = "/") #create proportion for each taxon/OTU within each sample
  sweep(temp, MARGIN = 1, STATS = rowMax(temp), FUN = "/")
}
#df <- as.data.frame(matrix(sample(1:10, 100, replace = T), ncol = 10, nrow = 10))
#eDNAINDEX(df)
