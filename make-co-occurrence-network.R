# A function to calculate a co-occurrence network from a phyloseq object

# Usage:
# ps : your phyloseq object
# corrlevel = 0.6
# alpha = 0.05
# make_co_occurence_network(ps, 0.3, 0.05)

library(phyloseq)
library(biomformat)
library(psych)
library(igraph)

make_co_occurrence_network <- function(phyloseq, corrlevel, alpha) {
  
  dat = data.frame(otu_table(phyloseq))
  
  # calculate correlations
  cor = corr.test(t(dat),use="pairwise",method="spearman",adjust="BH",alpha=alpha)
  cor.r= cor$r 
  cor.p= cor$p 
  cor.r[cor.p>alpha|cor.r<corrlevel] = 0 # select rho and pval
  # remove non-significant OTUs (e.g., were assigned "0" in prev step)w
#  cor.r = cor.r[colSums(cor.r) > 0, rowSums(cor.r) > 0]
  # convert the adjacency matrix to a sparse matrix for a faster and more memory-efficient analysis
  cor.r = Matrix(cor.r, sparse = T)
  
  # make network
  net = graph_from_adjacency_matrix(cor.r,mode="undirected", diag=F)
  
  # # weights cannot be negative for clustering
  # net.weight = E(net)$weight
  # E(net)$weight = NA
  # E(net)$width = abs(net.weight)
  
  #  modularity
  fc = cluster_fast_greedy(net,weights =NULL)
  modularity = modularity(net,membership(fc))
  comps = membership(fc)
  
  # random network
  random.net =erdos.renyi.game(length(V(net)),length(E(net)),type="gnm")
  
  # graph properties
  net.charac=c( transitivity(net), # clustering.coefficient
                transitivity(random.net), # clustering coefficient of random networks
                modularity,# modularity
                mean(igraph::degree(net)),# mean degree
                length(E(net)),# size - number of edges
                length(V(net)),# order - number of vertices
                edge_density(net,loops=FALSE),#edge.density
                mean_distance(net),#average path length
                no.clusters(net),#number.cluster
                mean(degree(net)/length(V(net))),#normalized degree
                centralization.betweenness(net)$centralization, #betweenness centrality
                average.path.length(net) # average shortest path
  )
  results = data.frame(property=c("clustering.coefficient", 
             "random.clustering.coefficient", 
             'modularity', 
             'mean.degree', 
             'size', 
             'order', 
             'edge.density', 
             'mean.distance', 
             'no.clusters',
             'norm.degree',
             'betweenness.centrality',
             'mean.shortest.path'),
             value=round(net.charac, 4)
             )
  return(results)
}
