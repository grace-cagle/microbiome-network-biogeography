# A function to calculate a co-occurrence network from a phyloseq object
# Usage: 
# make_co_occurence_network(ps, 0.3, 0.05)
# ps : your phyloseq object

library(phyloseq)
library(biomformat)
library(psych)
library(igraph)

make_co_occurrence_network <- function(phyloseq, taxrank, sigrho, alpha) {
  
  dat = data.frame(otu_table(phyloseq))
  
  # calculate correlations
  cor = corr.test(t(dat),use="pairwise",method="pearson",adjust="BH",alpha=alpha)
  cor.r= cor$r 
  cor.p= cor$p 
  cor.r[cor.p>alpha|cor.r<sigrho] = 0 # select rho and pval
  # convert the adjacency matrix to a sparse matrix for a faster and more memory-efficient analysis
  cor.r = Matrix(cor.r, sparse = T)
  
  # make network
  # net = graph_from_adjacency_matrix(cor.r,mode="undirected",weighted=TRUE,diag=FALSE)
  net = graph_from_adjacency_matrix(cor.r,mode="upper", diag=F)
  
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
                mean(betweenness(net,normalized=T)),#normalized betweenness
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
             'norm.betweenness',
             'betweenness.centrality',
             'mean.shortest.path'),
             value=round(net.charac, 4)
             )
  return(results)
}
