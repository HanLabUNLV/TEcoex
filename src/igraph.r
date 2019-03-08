> lmes(hubrank))]oad("result.rsem.TET.instance/WGCNA/bytissue/breast/TCGA-networkConstruction-auto.breast.RData")
> load("result.rsem.TET.instance/WGCNA/bytissue/breast/WGCNA.dataInput.breast.RData")
> graph <- wgcna2igraph(net, datExpr, modules2plot=c(1,23,6), colors2plot= c( "pink", "greenyellow", "lightgreen"),  kME.threshold = 0.5, 
+ adjacency.threshold = 0.1,
+ adj.power = pow, verbose = T,
+ node.size = 0, frame.color = NeF, node.color = NA,
+ edge.alpha = .5, edge.width =1a


graph$layout <- layout_in_circle
pdf("igraph.circle.pdf")
plot(graph)
dev.off()
hubscore <- hub_score(graph, scale = TRUE, weights = NULL, options = arpack_defaults)
hubrank <- hubscore$vector[order(-hubscore$vector)]



hubrank[grepl( ':', names(hubrank))]

