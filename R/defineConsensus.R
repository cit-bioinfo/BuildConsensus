# Rd
# description >> Compare and cluster various classification systems using a "Cluster of Cluster" approach.
# argument
# item >> annot.class >> Dataframe of samples annotated according to the several classification systems to compare.
# item >> outdir >> Path to the directory where to store plots and results. A new directory will be created if the supplied path does not exist.
# item >> maxK >> Integer value. maximum cluster number to evaluate.
# value >> A list of length maxK. Each element is a list containing consensusMatrix (numerical matrix), consensusTree (hclust), consensusClass (consensus class asssignments)
# author >> Aurelie Kamoun
# keyword >> methods
# details >> ...
# seealso >> ...
# references >> ...
# examples >> ...
# end
consensus.CoC <- function(annot.class, outdir, maxK = 8, reps = 100, seed = 12,
                          innerLinkage = "ward.D2", finalLinkage = "ward.D2", distance = "euclidean", ...){
  require(ConsensusClusterPlus)
  require(ComplexHeatmap)
  
  cat("You can have a 2 minutes break \n")
  
  if(!dir.exists(outdir)) {dir.create(outdir)}
  
  data <- apply(annot.class, 2, function(classif){sapply(unique(classif), function(cl){classif == cl})})
  cluster.mat <- t(do.call(cbind, data)) 
  
  rownames(cluster.mat) <- paste(rep(sapply(colnames(annot.class), function(x){unlist(strsplit(x, split = "\\."))[[1]]}), 
                                     times = apply(annot.class, 2, function(cl){length(table(cl))})),
                                 unlist(apply(annot.class, 2, function(cl) unique(cl))), sep = "_")
  
  coc <- ConsensusClusterPlus(cluster.mat, maxK = maxK, reps = reps,
                              title = outdir, distance = distance,
                              innerLinkage = innerLinkage, finalLinkage = finalLinkage,
                              corUse = "complete.obs",
                              seed = seed, plot = "pdf", ...)
  
  save(coc, file = file.path(outdir, "coc_results.RData"))

  icl <- calcICL(coc, title = outdir, plot = "pdf")
  
  pdf(file = file.path(outdir, "Heatmaps.pdf"), width = 12, height = 6)
  for(k in 2:maxK){
    col.class <- coc[[k]]$clrs[[3]]
    names(col.class) <- c(1:k)
    annotSamples <- HeatmapAnnotation(data.frame(consensus = coc[[k]]$consensusClass),
                                      col = list(consensus = col.class)
    )
    ComplexHeatmap::draw(Heatmap(1*cluster.mat, col = c("ghostwhite", "dimgrey"), name = "",
                  row_names_gp = gpar(fontsize = 10), 
                  show_column_names = F, top_annotation = annotSamples,
                  cluster_columns = coc[[k]]$consensusTree,
                  clustering_method_rows = "ward.D2")
    )
  }
  dev.off() 
  
  # plot average silhouette width
  sil.mean <- sapply(coc[2:length(coc)], function(k){mean(summary(silhouette(k$consensusClass,1-k$consensusMatrix))$clus.avg.widths)})
  names(sil.mean) <- 2:length(coc)
  
  point.col <- c("black", "red")[as.numeric(as.factor(sil.mean == max(sil.mean)))]
  
  pdf(file = file.path(outdir, "silhouette_vs_clusterNb.pdf"), height = 5, width = 5)
  plot(as.numeric(names(sil.mean)), sil.mean, col = point.col, pch = 18, ylim = c(min(sil.mean)-.01,max(sil.mean)+.01),
       xlab = "Nb of clusters", ylab = "Mean Silhouette Width")
  
  dev.off()
  
  invisible(coc)
}


# Rd
# description >> Compare and cluster various classification systems applying the same method as in Guinney et al (see reference).
# argument
# item >> annot.class >> Dataframe of samples annotated according to the several classification systems to compare.
# item >> I.values >> A vector of inflation factors values to use with MCL algorithm. The function running time linearly depends on the given number of values (about 3 minutes for each inflation factor value)
# item >> outdir >> Path to the directory where to store plots and results. A new directory will be created if the supplied path does not exist.
# item >> resamp >> Value between 0 and 1. Proportion of samples to use for MCL bootstrap iterations. Default is 0.8.
# item >> n.iter >> Number of MCL bootstrap iterations for each inflation factor value. Default is 1000.
# item >> pval.cut >> Cut-off for adjusted P-value to select significant network edges (Fisher test). Default is 0.001 
# value >> A list with the same length of I.values numerical vector. Each element is a list containing the cluster assignments (cl), the co-classification matrix (consensusMat), a subtype stability measure (subtype.stab), the mean weighted average silhouette width (wsil.mean)
# author >> Aurelie Kamoun
# keyword >> methods
# details >> ...
# seealso >> ...
# references >> 1.Guinney, J. et al. The consensus molecular subtypes of colorectal cancer. Nat Med advance online publication, (2015).
# examples >> ...
# end
consensus.MCL <- function(annot.class, I.values = 3.2, outdir, resamp = .8, n.iter = 1000, pval.cut = 0.001, seed = 42,
                          sim.method = c("Jaccard", "CohenKappa")[1], filter.ini = c("fisher", "cohenkappa")[1], ck.cut = .2){
  
  require(proxy)
  require(cluster)
  
  cat("Final results in ",  round(.4 * n.iter * length(I.values) / 60) , " minutes", "\n")
  
  mclfilesdir <- file.path(outdir ,"MCL_files")
  
  if(!dir.exists(outdir)) {dir.create(outdir)}
  if(!dir.exists(mclfilesdir)) {dir.create(mclfilesdir)}
  
  data <- apply(annot.class, 2, function(classif){sapply(unique(classif), function(cl){classif == cl})})
  cluster.mat <- t(do.call(cbind, data)) 
  
  rownames(cluster.mat) <- paste(rep(sapply(colnames(annot.class), function(x){unlist(strsplit(x, split = "\\."))[[1]]}), 
                                     times = apply(annot.class, 2, function(cl){length(table(cl))})),
                                 unlist(apply(annot.class, 2, function(cl) unique(cl))), sep = "_")
 

# Build network based on Jaccard similarities and prepare data for MCL
  
  if(sim.method == "Jaccard"){
    # Compute Jaccard similarity matrix
    S <- as.matrix(1 - dist(cluster.mat, method = "Jaccard"))
    S[lower.tri(S)] <- 0
    E <- which(S > 0, arr.ind = T, useNames = T)
  } else if(sim.method == "CohenKappa"){
    # Compute Cohen's kappa similarity matrix
    CohenKappa. <- function(i, j, data) {CohenKappa(data[, i], data[, j])}
    CohenKappa.v <- Vectorize(CohenKappa., vectorize.args  = list("i", "j"))
    S <- outer(1:nrow(cluster.mat), 1:nrow(cluster.mat), CohenKappa.v, data = t(cluster.mat))
    dimnames(S) <- list(rownames(cluster.mat), rownames(cluster.mat))
    S[which(S == 1)] <- 0
    E <- which(1 * upper.tri(S) * S > 0, arr = T)
  }
  network <- data.frame(node1 = rownames(S)[E[,1]], node2 = rownames(S)[E[,2]], weight = S[E])
  
  if(filter.ini == "fisher"){
    network$pval <- apply(network[, 1:2], 1, function(x){fisher.test(table(cluster.mat[as.character(x[1]), ],  
                                                                          cluster.mat[as.character(x[2]), ]), 
                                                                    alternative = "greater")$p.value})
    network$pval.adj <- p.adjust(network$pval, method = "BH")
    write.table(network[which(network$pval.adj < pval.cut), c("node1", "node2", "weight")], file = file.path(mclfilesdir, "data.abc"), sep = " ", row.names = F, col.names = F, quote = F)
  } else if(filter.ini == "cohenkappa"){
    if(sim.method == "CohenKappa"){
      write.table(network[which(network$weight >= ck.cut), c("node1", "node2", "weight")], file = file.path(mclfilesdir, "data.abc"), sep = " ", row.names = F, col.names = F, quote = F)
    } else{
      CohenKappa. <- function(i, j, data) {CohenKappa(data[, i], data[, j])}
      CohenKappa.v <- Vectorize(CohenKappa., vectorize.args  = list("i", "j"))
      S. <- outer(1:nrow(cluster.mat), 1:nrow(cluster.mat), CohenKappa.v, data = t(cluster.mat))
      dimnames(S.) <- list(rownames(cluster.mat), rownames(cluster.mat))
      S.[which(S. == 1)] <- 0
      E. <- which(1 * upper.tri(S.) * S. > 0, arr = T)
      ck.score <- S.[E.]
      write.table(network[which(ck.score >= ck.cut), c("node1", "node2", "weight")], file = file.path(mclfilesdir, "data.abc"), sep = " ", row.names = F, col.names = F, quote = F)
    }
  }
      
  system(paste("cd ", mclfilesdir, "; mcxload -abc data.abc --stream-mirror -write-tab data.tab -o data.mci", sep = ""), ignore.stdout = T, ignore.stderr = T)

      
# Partition network using MCL algorithm
    
  MCL <- list()
    
  for(I in I.values){            
    
    system(paste("cd ", mclfilesdir, "; mcl data.mci -I ", I, sep = ""), ignore.stdout = T, ignore.stderr = T)
    system(paste("cd ", mclfilesdir, "; mcxdump -icl out.data.mci.I", I * 10," -tabr data.tab -o dump.data.mci.I", I * 10, sep = ""), ignore.stdout = T, ignore.stderr = T)
    
    res.all <- strsplit(scan(file = file.path(mclfilesdir, paste("dump.data.mci.I", 10 * I, sep = "")), what = "", sep = "\n", quiet = T), "\t")
    cl.all <- rep(1:length(res.all), times = lapply(res.all, length))
    names(cl.all) <- unlist(res.all)
    
    MCL[[paste(I)]]$cl <- cl.all
    
    network.all <- crossprod(table(cl.all[rownames(cluster.mat)],rownames(cluster.mat)))
    network.all <- network.all[rownames(cluster.mat), rownames(cluster.mat)]
    

# Initialize co-classification matrix and subtypes stability measures 
    
    consensusMat <- matrix(0, ncol = nrow(cluster.mat), nrow = nrow(cluster.mat),
                           dimnames = list(rownames(cluster.mat), rownames(cluster.mat)))
    
    network.stab <- matrix(0, ncol = nrow(cluster.mat), nrow = nrow(cluster.mat),
                           dimnames = list(rownames(cluster.mat), rownames(cluster.mat)))

    
# Bootstrap iterations: iterate network contruction and partitioning using only a subset of samples
    
    if(!dir.exists(file.path(mclfilesdir, "tmp"))) {dir.create(file.path(mclfilesdir, "tmp"))}
    
    cat("Iterating for I = ", I, "\n")
    
    set.seed(seed)
    for(i in 1:n.iter){
      
      samp <- sample(1:ncol(cluster.mat), round(resamp * ncol(cluster.mat)))
      dat <- cluster.mat[, samp]
      
      if(sim.method == "Jaccard"){
        S <- as.matrix(1 - dist(dat, method = "Jaccard"))
        S[lower.tri(S)] <- 0
        E <- which(S > 0, arr.ind = T, useNames = T)
      } else if(sim.method == "CohenKappa"){
        S <- outer(1:nrow(dat), 1:nrow(dat), CohenKappa.v, data = t(dat))
        dimnames(S) <- list(rownames(dat), rownames(dat))
        S[which(S == 1)] <- 0
        E <- which(1 * upper.tri(S) * S > 0, arr = T)
      }
      
      network <- data.frame(node1 = rownames(S)[E[,1]], node2 = rownames(S)[E[,2]], weight = S[E])
      
      if(filter.ini == "fisher"){
        network$pval <- apply(network[, 1:2], 1, function(x){fisher.test(table(dat[as.character(x[1]), ],  
                                                                            dat[as.character(x[2]), ]), 
                                                                      alternative = "greater")$p.value})
        network$pval.adj <- p.adjust(network$pval, method = "BH")
        write.table(network[which(network$pval.adj < pval.cut), c("node1", "node2", "weight")], file = file.path(mclfilesdir, "tmp/iter.abc"), sep = " ", row.names = F, col.names = F, quote = F)
      } else if(filter.ini == "cohenkappa"){
        if (sim.method == "CohenKappa"){
          write.table(network[which(network$weight >= ck.cut), c("node1", "node2", "weight")], file = file.path(mclfilesdir, "tmp/iter.abc"), sep = " ", row.names = F, col.names = F, quote = F)
        } else {
          CohenKappa. <- function(i, j, data) {CohenKappa(data[, i], data[, j])}
          CohenKappa.v <- Vectorize(CohenKappa., vectorize.args  = list("i", "j"))
          S. <- outer(1:nrow(dat), 1:nrow(dat), CohenKappa.v, data = t(dat))
          dimnames(S.) <- list(rownames(dat), rownames(dat))
          S.[which(S. == 1)] <- 0
          E. <- which(1 * upper.tri(S.) * S. > 0, arr = T)
          ck.score <- S.[E.]
          write.table(network[which(ck.score >= ck.cut), c("node1", "node2", "weight")], file = file.path(mclfilesdir, "tmp/iter.abc"), sep = " ", row.names = F, col.names = F, quote = F)
        }
      }
      
      system(paste("cd ", file.path(mclfilesdir, "tmp"), "; mcxload -abc iter.abc --stream-mirror -write-tab iter.tab -o iter.mci", sep = ""), ignore.stdout = T, ignore.stderr = T)
      system(paste("cd ", file.path(mclfilesdir, "tmp"), "; mcl iter.mci -I ",I, sep = ""), ignore.stdout = T, ignore.stderr = T)
      system(paste("cd ", file.path(mclfilesdir, "tmp"), "; mcxdump -icl out.iter.mci.I", I * 10," -tabr iter.tab -o dump.iter.mci.I", I * 10, sep = ""), ignore.stdout = T, ignore.stderr = T)
      
      res <- strsplit(scan(file = file.path(mclfilesdir, paste("tmp/dump.iter.mci.I", 10 * I, sep = "")), what = "", sep = "\n", quiet = T), "\t")
      cl <- rep(1:length(res), times = lapply(res, length))
      names(cl) <- unlist(res)
      
      assocMat <- crossprod(table(cl[rownames(consensusMat)], rownames(consensusMat)))
      assocMat <- assocMat[rownames(consensusMat), rownames(consensusMat)]
      consensusMat <- consensusMat + assocMat
      
      network.stab <- network.stab + apply(assocMat == network.all, 1, as.numeric)
      
    }

    
# Bootstrap aggregating 
    
    MCL[[paste(I)]]$consensusMat <- consensusMat/n.iter
    MCL[[paste(I)]]$subtype.stab <- apply(network.stab, 2, mean)/n.iter
    
    x <- rep(1:length(res.all), times = unlist(lapply(res.all, length)))
    D <- 1 - MCL[[paste(I)]]$consensusMat[unlist(res.all), unlist(res.all)]
    
    sil <- silhouette(x, D)
    sil.val <- sil[, 3]
    names(sil.val) <- unlist(res.all)
    
    w <- MCL[[paste(I)]]$subtype.stab[unlist(res.all)]
    
    MCL[[paste(I)]]$wsil.mean <- mean(sil.val*w)
    
  } # repeat for each inflation factor value in I.values
  
  cat("Plotting results \n")
  
  
  #Silhouette plot vs inflation factor
  
  wsil.mean <- unlist(lapply(MCL, function(x){x$wsil.mean}), use.names = T)
  point.col <- c("black", "red")[as.numeric(as.factor(wsil.mean == max(wsil.mean)))]
  point.cl <- unlist(lapply(MCL, function(x){length(unique(x$cl))}), use.names = T)
  
  pdf(file = file.path(outdir, "silhouette_vs_InflationFactor.pdf"), height = 5, width = 5)
  plot(as.numeric(names(wsil.mean)), wsil.mean, col = point.col, pch = 18, ylim = c(min(wsil.mean) - .01,max(wsil.mean) + .01),
       xlab = "Inflation factor", ylab = "Mean Weighted Silhouette Width")
  text(as.numeric(names(wsil.mean)), wsil.mean, labels = point.cl, pos = 3, cex = .7)
  dev.off()
  
  I.opt <- names(which(wsil.mean == max(wsil.mean)))
  write.table(as.data.frame(sapply(I.opt, function(x){MCL[[x]]$cl})), 
              file = file.path(outdir, "optimalI_results.txt"),row.names = T, quote = F, sep = "\t")
  
  #Heatmaps plot
  
  pdf(file = file.path(outdir, paste("Heatmap_I", 10 * as.numeric(I), ".pdf", sep = "")), width = 5, height = 5)
  layout(matrix(1:6, ncol = 2, byrow = T))
  
  for (I in names(MCL)){
    ComplexHeatmap::draw(Heatmap(MCL[[I]]$consensusMat, col = c("white", "orange"), show_heatmap_legend = F, 
                                 row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
                                 column_title = paste("I = ", I, sep = "")))
  }
  
  dev.off()
  
  #Iwrite <- if(length(I.values) > 1) paste(10  *I.values[1], "-", 10 * I.values[length(I.values)], sep = "") else as.character(10 * I.values[1])
  save(MCL, file = file.path(outdir, "mcl_results.RData"))
       
  invisible(MCL)
  

}


# Rd
# description >> Compare and cluster various classification systems using Cohen's Kappa.
# argument
# item >> annot.class >> Dataframe of samples annotated according to the several classification systems to compare.
# item >> outdir >> Path to the directory where to store plots and results. A new directory will be created if the supplied path does not exist.
# item >> CohenKappa.cut >> Value of Cohen's Kappa to select significant associations.
# value >> A list of length 3 containing the Cohen's Kappa measures for any pair of subtypes (CohenKappaMat), the corresponding network generated (graph), the cluster assignments (cl)
# author >> Aurelie Kamoun
# keyword >> methods
# details >> ...
# seealso >> ...
# references >> ...
# examples >> ...
# end
consensus.CIT <- function(annot.class, outdir, CohenKappa.cut = seq(0.1, .7, .1)){
  require(cluster)
  require(DescTools) 
  require(igraph)
  
  if(!dir.exists(outdir)) {dir.create(outdir)}
  
  # Binarize subtypes data
  data <- apply(annot.class, 2, function(classif){sapply(unique(classif), function(cl){classif == cl})})
  cluster.mat <- t(do.call(cbind, data)) 
  
  rownames(cluster.mat) <- paste(rep(sapply(colnames(annot.class), function(x){unlist(strsplit(x, split = "\\."))[[1]]}), 
                                     times = apply(annot.class, 2, function(cl){length(table(cl))})),
                                 unlist(apply(annot.class, 2, function(cl) unique(cl))), sep = "_")
  
  # Compute similarity matrix (Cohen Kappa)
  CohenKappa. <- function(i, j, data) {CohenKappa(data[, i], data[, j])}
  CohenKappa.v <- Vectorize(CohenKappa., vectorize.args  = list("i", "j"))
  S <- outer(1:nrow(cluster.mat), 1:nrow(cluster.mat), CohenKappa.v, data = t(cluster.mat))
  dimnames(S) <- list(rownames(cluster.mat), rownames(cluster.mat))
  S[which(S == 1)] <- 0
  
  CIT <- list()
  CIT$CohenKappaMat <- S
  
  mod <- numeric()
  ncl <- numeric()
  
  for(ck in CohenKappa.cut){
    
    # Resulting network
    E <- which(1*upper.tri(S)*S > ck, arr = T)
    network <- data.frame(node1 = rownames(S)[E[,1]], node2 = rownames(S)[E[,2]], weight = S[E])
  
    # Visualize graph and optimal community structure (maximizing modularity)
    CIT[[paste(ck)]]$graph <- graph.data.frame(network, directed=F)
    CIT[[paste(ck)]]$cl <- membership(cluster_optimal(CIT[[paste(ck)]]$graph))
  
    color.nodes <- rainbow(length(unique(CIT[[paste(ck)]]$cl)))[CIT[[paste(ck)]]$cl]
    names(color.nodes) <- names(CIT[[paste(ck)]]$cl)
  
    V(CIT[[paste(ck)]]$graph)$color <- color.nodes[names(V(CIT[[paste(ck)]]$graph))]
    V(CIT[[paste(ck)]]$graph)$shape <- "square"
  
    mod <- append(mod, modularity(cluster_optimal(CIT[[paste(ck)]]$graph)))
    ncl <- append(ncl, length(unique(CIT[[paste(ck)]]$cl)))
    
    pdf(file = file.path(outdir, paste("plot_graph_CKcut", ck, ".pdf")), width = 6, height = 5)
    par(mar = c(4, 3, 3, 10), cex.main = .9)
    plot(CIT[[paste(ck)]]$graph, layout = layout.graphopt, edge.width = 10 * network$weight,
         main = paste("modularity : ", modularity(cluster_optimal(CIT[[paste(ck)]]$graph)), sep = ""))
    legend("bottomright", legend = seq(.2, .8, .2), col = "gray65", lwd = seq(2, 8, 2), bty = "n", title = "Cohen's kappa", inset = c(-.5, 0), cex = .8)
    dev.off()
  }
  
  pdf(file = file.path(outdir, "Modularity_vs_CohenKappa.pdf"), height = 5, width = 5)
  point.col <- c("black", "red")[as.numeric(as.factor(mod == max(mod)))]
  plot(CohenKappa.cut, mod, col = point.col, pch = 18, xlab = "Cohen-Kappa cut-off", ylab = "Optimal graph modularity",
       ylim = c(0, max(mod) + .1))
  text(CohenKappa.cut, mod, labels = ncl, pos = 3, cex = .9)
  dev.off()
  
  save(CIT, file = file.path(outdir, paste("cit_results_CKcut", CohenKappa.cut[1], "-", ck, ".RData", sep = "")))
  
  invisible(CIT)

}

#consensus.compare <- function(mcl_results = NULL, coc_results = NULL, cit_results = NULL){}
