##### Function to calculate eDNA index
eDNAINDEX <- function(x) { #where x is a dataframe with taxa/OTUs/etc in rows, and samples in columns
  rowMax <- function(x){apply(x, MARGIN = 1, FUN = max)}
  temp <- sweep(x, MARGIN = 2, STATS = colSums(x), FUN = "/") #create proportion for each taxon/OTU within each sample
  sweep(temp, MARGIN = 1, STATS = rowMax(temp), FUN = "/")
}
#df <- as.data.frame(matrix(sample(1:10, 100, replace = T), ncol = 10, nrow = 10))
#eDNAINDEX(df)
