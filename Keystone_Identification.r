######################################################################
### Calculate correlation coefficients between microbial abundances
######################################################################
library(Hmisc)

# Read the otu-sample matrix, the row is sample, and the column is otu
otu <- read.csv(file.choose(), head=T, row.names=1)

# Calculate correlation
sp.cor<-rcorr(t(otu),type="spearman")

# Extract r and p value matrix
occor.r <- sp.cor$r
occor.p <- sp.cor$p

# Use the Benjamini-Hochberg ("FDR-BH") method for multiple testing correction
p <- p.adjust(occor.p, method="BH")

# Determine the threshold for the existence of interaction relationships between species, and convert non-compliant data in the correlation R matrix to 0
occor.r[occor.p>0.01|abs(occor.r)<0.8] = 0
diag(occor.r) <- 0

# Save occor.r as csv file
# write.csv(occor.r,file="Correlation calculation results.csv")

# Keep data based on r-values and p-values filtered above
z <- occor.r * occor.p
diag(z) <- 0 # Convert the values in the diagonal of the correlation matrix (representing autocorrelation) to 0
head(z)[1:6,1:6]
z[abs(z)>0]=1
z
adjacency_unweight <- z


######################################################################
### igraph package computing network module
######################################################################
library(igraph)
# Input data example, adjacency matrix
# This is a microbial interaction network. The value "1" indicates that there is interaction between microbial OTUs, and "0" indicates that there is no interaction
head(adjacency_unweight)[1:6] # Adjacency matrix type network file

# Adjacency matrix -> igraph adjacency list to obtain a non-weighted undirected network
igraph <- graph_from_adjacency_matrix(as.matrix(adjacency_unweight), mode = 'undirected', weighted = NULL, diag = FALSE)
igraph # Adjacency list for igraph

# Calculate node degree
V(igraph)$degree <- degree(igraph)

# Module division, details cluster_fast_greedy, there are multiple models
set.seed(123)
V(igraph)$modularity <- membership(cluster_fast_greedy(igraph))

# Output the list of each node name, node degree, and the modules it is divided into
nodes_list <- data.frame(
  nodes_id = V(igraph)$name, 
  degree = V(igraph)$degree, 
  modularity = V(igraph)$modularity
)
head(nodes_list) # Node list, including node name, node degree, and the modules it is divided into

# write.table(nodes_list, 'nodes_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)



######################################################################
### Calculate intra-module connectivity (Zi) and inter-module connectivity (Pi)
######################################################################
source('Zi-Pi_calculation.r')

# The above-mentioned adjacency matrix type network file
adjacency_unweight 

# Node attribute list, including the modules divided by the node
nodes_list <- read.delim('nodes_list.txt', row.names = 1, sep = '\t', check.names = FALSE)

# The node order of the two files must be consistent
nodes_list <- nodes_list[rownames(adjacency_unweight), ]

# Calculate intra-module connectivity (Zi) and inter-module connectivity (Pi)
# Specify the column names of the adjacency matrix, node list, node degree and module degree in the node list
zi_pi <- zi.pi(nodes_list, adjacency_unweight, degree = 'degree', modularity_class = 'modularity')
head(zi_pi)

# write.table(zi_pi, 'zi_pi_result.txt', sep = '\t', row.names = FALSE, quote = FALSE)

######################################################################
### Nodes can be divided into 4 types according to the threshold, and a graph can be drawn to show their distribution.
######################################################################
library(ggplot2)

zi_pi <- na.omit(zi_pi)
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'
# write.csv(zi_pi,"zipi_results.csv")
ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.8, size = 6,shape=17) +
  scale_y_continuous(limits=c(-2,3)) +
  scale_color_manual(values = c("#8491B4FF","#91D1C2FF","#F39B7FFF", "#4DBBD5FF"), 
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs')) +
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62,linetype=2,size=1) +
  geom_hline(yintercept = 2.5,linetype=2,size=1) +  
  theme_bw() +
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")))
