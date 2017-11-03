## Heteromotility plotting and analysis functions
source('sigbra.R')
library(devtools)
#source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

## Plotting

pca_plots <- function(df, df.class, df.plate, plot_path, class_num, clust_pcs=10, prefix=NULL){
  # Generates 2D PCA plots with a defined set of labels and k-means clusters.
  # 
  # Parameters
  # ----------
  # df : data frame, N x M. samples as rows, features as columns.
  # df.class : vector numeric. N x 1. class labels.
  # df.plate : vector numeric. N x 1. plate number labels.
  # plot_path : string. directory for plot outputs.
  # class_num : integer. k argument for k-means clustering.
  # clust_pcs : integer. number of PCs to use for k-means clustering.
  # prefix : character vector. prefix for filenames. Default = None.
  # 
  # Returns
  # -------
  # df.comp : data frame. N x 30 of principal component scores.
  # 
  require(ggplot2)
  # PCA
  df.pca_eigens <- prcomp(df)
  df.comp <- data.frame(df.pca_eigens$x[,1:30])
  df.k <- kmeans(df.comp[,1:clust_pcs], class_num, nstart=25, iter.max=1000)
  
  # plot PC1 vs PC2
  pca_theme <- theme(axis.text.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks=element_blank(),
                     legend.position="right",
                     panel.background=element_blank(),
                     panel.border=element_blank(),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),
                     plot.background=element_blank())
  pca_k <- ggplot(df.comp, aes(x = PC1, y = PC2, col = as.factor(df.k$cluster))) + geom_point(size=1) + pca_theme + scale_colour_discrete(name = "Cluster")
  ggsave(pca_k, filename = paste(prefix, "pca_k.png", sep=''), path = plot_path, width = 4, height = 4, units = "in")
  pca_class <- ggplot(df.comp, aes(x = PC1, y = PC2, col = as.factor(df.class))) + geom_point(size=1) + pca_theme + scale_colour_discrete(name = "Class")
  ggsave(pca_class, filename = paste(prefix, "pca_class.png", sep=''), path = plot_path, width = 4, height = 4, units = "in")
  pca_plate <- ggplot(df.comp, aes(x = PC1, y = PC2, col = as.factor(df.plate))) + geom_point(size=1) + pca_theme + scale_colour_discrete(name = "Plate")
  ggsave(pca_plate, filename = paste(prefix, "pca_plate.png", sep=''), path = plot_path, width = 4, height = 4, units = "in")
  return(df.comp)
}

tsne_plots <- function(df, df.groups, df.class, df.plate, plot_path, df.comp = NULL, class_num = NULL, experiment = NULL, method='ward.D2', perplexity=NULL, max_iter=5000, width=4, height=4){
  # Generates tSNE plots with cluster, class, and plate labels.
  # 
  # Parameters
  # ----------
  # df : data frame, N x M. samples as rows, features as columns.
  # df.groups : vector numeric. N x 1. cluster labels.
  # df.class : vector numeric. N x 1. class labels.
  # df.plate : vector numeric. N x 1. plate number labels.
  # plot_path : string. directory for plot outputs.
  # df.comp : data frame. N x D matrix of principal component scores, if plotting
  #           of both raw and PCA values is desired.
  # class_num : integer. k parameter for h-clustering of PCA values. 
  # experiment : string. Used as the prefix title for plot exports.
  # method : string. hierarchical clustering method for PCA values.
  # perplexity : vector numeric. perplexity values to try. 
  #              if NULL, a default set [10,...,70] is used.
  # max_iter : integer. maximum iterations of the tSNE Barnes-Hut algorithm.
  # width : integer. inches width of plots.
  # height : integer. inches height of plots.
  #
  # Returns
  # -------
  # None.
  #
  require(Rtsne)
  require(ggplot2)
  no_axes_bg <- theme(axis.line=element_blank(),
                      axis.text.x=element_blank(),
                      axis.text.y=element_blank(),
                      axis.ticks=element_blank(),
                      axis.title.x=element_blank(),
                      axis.title.y=element_blank(),
                      legend.position="none",
                      panel.background=element_blank(),
                      panel.border=element_blank(),
                      panel.grid.major=element_blank(),
                      panel.grid.minor=element_blank(),
                      plot.background=element_blank())
  
  no_axes_w_legend <- theme(axis.line=element_blank(),
                            axis.text.x=element_blank(),
                            axis.text.y=element_blank(),
                            axis.ticks=element_blank(),
                            axis.title.x=element_blank(),
                            axis.title.y=element_blank(),
                            legend.position="right",
                            panel.background=element_blank(),
                            panel.border=element_blank(),
                            panel.grid.major=element_blank(),
                            panel.grid.minor=element_blank(),
                            plot.background=element_blank())
  
  if (is.null(perplexity)){
    pset = c(10, 20, 30, 50, 70)
  } else { 
    pset = c(perplexity) 
    }
  
  costs <- numeric(length(pset))
  
  i <- 1
  for (perplexity in pset){
    df.tsne <- Rtsne(df, perplexity = perplexity, max_iter=max_iter)
    tsne_class <- ggplot(data.frame(df.tsne$Y), aes(x=X1, y=X2, col = as.factor(df.class))) + geom_point(size=1) + labs(title = paste(experiment, 'tSNE', sep = ' ')) + no_axes_bg
    ggsave(tsne_class, path = plot_path, filename = paste(experiment,"tsne_class_p", perplexity, ".png", sep=''), width = width, height = height, units = "in")
    tsne_wardD <- ggplot(data.frame(df.tsne$Y), aes(x=X1, y=X2, col = as.factor(df.groups))) + geom_point(size=1) + labs(title = paste(experiment, 'tSNE', sep = ' ')) + no_axes_bg
    ggsave(tsne_wardD, path = plot_path, filename = paste(experiment,"tsne_hclust_p", perplexity, ".png", sep=''), width = width, height = height, units = "in")
    tsne_plate <- ggplot(data.frame(df.tsne$Y), aes(x=X1, y=X2, col = as.factor(df.plate))) + geom_point() + labs(title = paste(experiment, 'tSNE', sep = ' ')) + no_axes_w_legend
    ggsave(tsne_plate, path = plot_path, filename = paste(experiment,"tsne_plate_p", perplexity, ".png", sep=''), width = width+3, height = height, units = "in") # extra width for legend
    
    costs[i] <- sum(df.tsne$costs)
  }
  
  if (!is.null(df.comp)){
    df.comp.tsne <- Rtsne(df.comp, perplexity = 50)
    
    df.comp.fit <- hclust(dist(df.comp), method = method)
    df.comp.groups <- cutree(df.comp.fit, k = class_num)
    tsne_pca_class <- ggplot(data.frame(df.comp.tsne$Y), aes(x=X1, y=X2, col = as.factor(df.class))) + geom_point(size=1) + labs(title = paste(experiment, 'tSNE', sep = ' ')) + no_axes_bg
    ggsave(tsne_pca_class, path = plot_path, filename = "tsne_pca_class.png", width = 4, height = 4, units = "in")
    tsne_pca_wardD <- ggplot(data.frame(df.comp.tsne$Y), aes(x=X1, y=X2, col = as.factor(df.comp.groups))) + geom_point(size=1) + labs(title = paste(experiment, 'tSNE', sep = ' ')) + no_axes_bg
    ggsave(tsne_pca_wardD, path = plot_path, filename = "tsne_pca_hclust.png", width = 4, height = 4, units = "in")
  }
}

## Plotting Colors and Themes

fgf_class_colors <- function(class_vector) {
  return_vector <- numeric(0)
  for (i in class_vector) {
    if (i == "FGF2+") {
      return_vector <- append(return_vector, 'blue')
    } else {
      return_vector <- append(return_vector, 'red')
    }
  }
  return(return_vector)
}

mef_class_colors <- function(class_vector) {
  return_vector <- numeric(0)
  for (i in class_vector) {
    if (i == "MycRas") {
      return_vector <- append(return_vector, 'purple')
    } else {
      return_vector <- append(return_vector, 'blue')
    }
  }
  return(return_vector)
}

## Clustering

som_plots <- function(df, grid_size, class_num, plot_path, method="ward.D2", iter=500){
  require(kohonen)
  som_grid <- somgrid(xdim = grid_size, ydim = grid_size, topo = "hexagonal") # set 20x20 hexagonal perceptron grid
  df.som_model <- som(as.matrix(df), grid=som_grid, rlen = iter, alpha = c(0.05, 0.01), keep.data = TRUE)
  require(RColorBrewer)
  somcol <- brewer.pal(11, 'RdBu')
  png(filename = paste(plot_path, "som_changes.png", sep = ''), width = 4, height = 4, units = "in", res = 600)
  plot(df.som_model, type = "changes")
  dev.off()
  
  png(filename = paste(plot_path, "som_counts.png", sep = ''), width = 4, height = 4, units = "in", res = 600)
  plot(df.som_model, type="count")
  dev.off()
  
  png(filename = paste(plot_path, "som_neighbor_distances.png", sep = ''), width = 4, height = 4, units = "in", res = 600)
  plot(df.som_model, type="dist.neighbours")
  dev.off()
  
  df.som_cluster <- cutree( hclust(dist(data.frame(df.som_model$codes)), method = method), k = class_num)
  pretty_palette <- c("#1f77b4", '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2')
  png(filename = paste(plot_path, "som_clusters.png", sep = ''), width = 4, height = 4, units = "in", res = 600)
  plot(df.som_model, type="mapping", bgcol = pretty_palette[df.som_cluster], main = "Clusters")
  add.cluster.boundaries(df.som_model, df.som_cluster)
  dev.off()
  
}

heatmap2_plots <- function(df, df.class, plot_path, class_num, method="ward.D2", prefix=''){
  # Heirarchical clustering
  require(gplots)
  df.d <- dist(df, method = "euclidean")
  df.fit <- hclust(df.d, method = method)
  df.groups <- cutree(df.fit, k = class_num)
  
  # make RowSideColors vector for heatmap.2
  rsc <- as.character(df.groups)
  
  hclust_fnx <- function(x) hclust(x, method=method)
  
  out_file <- paste(plot_path, prefix, "heatmap_groups_only.png", sep ='')
  png(out_file, width = 10, height = 10, units = "in", res = 300)
  heatmap.2(as.matrix(df), hclustfun = hclust_fnx, scale = "row", trace = 'none', col=redgreen, RowSideColors = rsc, margins = c(10,10))
  dev.off()
}


heatmap3_plots <- function(df, df.class, plot_path, class_num, method="ward.D2", prefix=''){
  # Heirarchical clustering
  require(gplots)
  df.d <- dist(df, method = "euclidean")
  df.fit <- hclust(df.d, method = method)
  df.groups <- cutree(df.fit, k = class_num)
  
  # make RowSideColors vector for heatmap.3
  r1 <- fgf_class_colors(df.class)
  r2 <- as.character(df.groups)
  rsc <- t( cbind(r1,r2) )
  rownames(rsc) <- c("Class", "Cluster")
  
  hclust_fnx <- function(x) hclust(x, method=method)
  
  out_file <- paste(plot_path, prefix, "heatmap_class.png", sep ='')
  png(out_file, width = 10, height = 10, units = "in", res = 300)
  heatmap.3(as.matrix(df), hclustfun = hclust_fnx, scale = "row", trace = 'none', col=redgreen, RowSideColors = rsc, margins = c(10,10))
  dev.off()
  
  return(df.groups)
}

## Cluster Membership

fgf_cluster_props <- function(class_df, df.groups){
  class1 = rep(0, max(df.groups))
  class2 = rep(0, max(df.groups))
  for (i in 1:nrow(class_df)){
    
    if (class_df$class[i] == 'FGF2+'){
      class1[df.groups[i]] <- class1[df.groups[i]] + 1
    }
    else {
      class2[df.groups[i]] <- class2[df.groups[i]] + 1
    }
  }
  
  props1 <- class1/sum(class1)
  props2 <- class2/sum(class2)
  
  return( list(fgf=class1, nofgf=class2, fgf_p = props1, nofgf_p = props2) )
}

fgf_plot_cluster_memberships <- function(cluster_membership, plot_path){
  d <- data.frame(rbind(cluster_membership$fgf_p, cluster_membership$nofgf_p))
  clust_nb <- ncol(d)
  names = c()
  for (i in 1:clust_nb){ names[i] = paste('Cluster', i, sep=' ') }
  colnames(d) <- names
  d$class <- c("FGF2+", "FGF2-")
  d.m <- melt(d, id.vars = "class")
  cluster_theme <- theme(
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank(),
    axis.line.x = element_line(size = 0.5),
    axis.line.y = element_line(size = 0.5),
    axis.text.x = element_text(size = 10), 
    text= element_text(size = 10)
  )
  clust_p <- ggplot(d.m, aes(variable, value, fill = class)) + geom_bar(position = "dodge", stat = "identity")
  clust_p <- clust_p + cluster_theme + scale_y_continuous(expand = c(0,0), limits = c(0, 1.10*max(d.m$value))) + labs(title = "Cluster Membership", y = "Proportion in Cluster", x = "")
  ggsave(clust_p, filename = "cluster_membership.png", path = plot_path, width = 4, height = 3, units = "in")
}

mef_cluster_props <- function(class_df, df.groups){
  class1 = rep(0, max(df.groups))
  class2 = rep(0, max(df.groups))
  for (i in 1:nrow(class_df)){
    
    if (class_df$class[i] == 'MycRas'){
      class1[df.groups[i]] <- class1[df.groups[i]] + 1
    }
    else {
      class2[df.groups[i]] <- class2[df.groups[i]] + 1
    }
  }
  
  props1 <- class1/sum(class1)
  props2 <- class2/sum(class2)
  
  return( list(mycras=class1, wt=class2, mycras_p = props1, wt_p = props2) )
}

mef_plot_cluster_memberships <- function(cluster_membership, plot_path, prefix=''){
  d <- data.frame(rbind(cluster_membership$mycras_p, cluster_membership$wt_p))
  names = c()
  for (i in 1:ncol(d)){ names[i] = paste('Cluster', i, sep=' ') }
  colnames(d) <- names
  d$class <- c("MycRas", "WT")
  d.m <- melt(d, id.vars = "class")
  cluster_theme <- theme(
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank(),
    axis.line.x = element_line(size = 0.5),
    axis.line.y = element_line(size = 0.5),
    axis.text.x = element_text(size = 10), 
    text= element_text(size = 10)
  )
  clust_p <- ggplot(d.m, aes(variable, value, fill = class)) + geom_bar(position = "dodge", stat = "identity")
  clust_p <- clust_p + cluster_theme + scale_y_continuous(expand = c(0,0), limits = c(0, 1.10*max(d.m$value))) + labs(title = "Cluster Membership", y = "Proportion in Cluster", x = "")
  ggsave(clust_p, filename = paste(prefix, "cluster_membership.png", sep=''), path = plot_path, width = 4, height = 3, units = "in")
}


density_cluster_membership <- function(df.groups, density.state, plot_path, prefix='', chi2test=T){
  n_groups <- max(df.groups)
  n_dens <- max(density.state)
  groups_mat <- matrix(nrow=n_groups, ncol=n_dens, 0)
  
  # Generate matrix of all proportions of each density state in each group
  for (i in 1:length(df.groups)){
    groups_mat[df.groups[i], density.state[i]] = groups_mat[df.groups[i], density.state[i]] + 1
  }
  
  density_contingency_table_preftest(groups_mat, plot_path, prefix=prefix)
  
  # Get proportion of density state in each cluster
  sums_mat <- t(replicate(n_groups, apply(groups_mat, 2, sum)))
  props_mat <- groups_mat/sums_mat
  
  props_df <- data.frame(props_mat)
  densclust_names <- c()
  for (i in 1:n_dens){densclust_names <- c(densclust_names, paste('Density ', i, sep=''))}
  colnames(props_df) <- densclust_names
  return(props_df)
}

plot_density_cluster_membership <- function(props_df, plot_path, prefix=''){
  # props_df : data.frame of State x Density_Cluster, containing proportion of each density cluster
  #           in each state.
  require(reshape2)
  cluster_theme <- theme(
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank(),
    axis.line.x = element_line(size = 0.5),
    axis.line.y = element_line(size = 0.5),
    axis.text.x = element_text(size = 10), 
    text= element_text(size = 10)
  )
  clust_names <- c()
  for (i in 1:nrow(props_df)){clust_names <- c(clust_names, paste('Cluster ', i, sep=''))}
  props_df$Cluster <- clust_names
  props.m <- melt(props_df, id.vars=c('Cluster'))
  clust_p <- ggplot(props.m, aes(variable, value, fill = as.factor(Cluster))) + geom_bar(position = "dodge", stat = "identity")
  clust_p <- clust_p + cluster_theme + scale_y_continuous(expand = c(0,0), limits = c(0, 1.10*max(props.m$value))) + labs(title = "Density ~ Cluster Membership", y = "Proportion in Cluster", x = "")
  ggsave(clust_p, path=plot_path, filename = paste(prefix, '_density_clusterprops.png', sep=''), height=3, width=4, units='in')
  
  clust_p2 <- ggplot(props.m, aes(Cluster, value, fill = variable)) + geom_bar(position = "dodge", stat = "identity")
  clust_p2 <- clust_p2 + cluster_theme + scale_y_continuous(expand = c(0,0), limits = c(0, 1.10*max(props.m$value))) + labs(title = "Cluster Membership ~ Density", y= 'Proportion in Cluster', x = "")
  ggsave(clust_p2, path=plot_path, filename = paste(prefix, '_cluster_densityprops.png', sep=''), height=3, width=4, units='in')
}

pca_loading_variation <- function(df, plot_path, prefix=''){
  # Saves PCA loadings and percentage variation captured at different levels of components
  df.pc <- princomp(df)
  aload <- abs( with(df.pc, unclass(loadings)) ) # get abs of all loadings
  pload <- sweep(aload, 2, colSums(aload), "/") # proportion per PC
  write.csv(x = pload, file = paste(plot_path, prefix, 'pca_loading_proportional.csv', sep=''), row.names = F)
  
  # get prop variance
  df.pca_eigens <- prcomp(df)
  s <- summary(df.pca_eigens)
  capture.output(s, file = paste(plot_path, prefix, 'pca_var_explained.txt', sep=''))
  write.csv(s$importance[3,], file = paste(plot_path, prefix, 'pca_cum_var_explained.csv', sep=''), row.names = F)
}

## Feature Specific Plots

add_sig_bars <- function(plot, plt_df, test_results, y_shift = 0.05, alpha = 0.05, text.size=5, lab.space=0.003){
  # Adds significance bars to a provided ggplot object
  #
  # Parameters
  # ----------
  # plot : ggplot geom_bar object.
  # plt_df : longform data frame with Mean and SE for a set of variables
  # test_results : wideform data frame with a column titled adj_p_val of p-values to use for sig bars
  # y_shift : float. Distance above highest error bar val to place sig bar.
  # alpha : float, (0, 1). significance level.
  # text.size : integer.
  # lab.space : float. Distance above bar to place sig mark.
  #
  # Returns
  # -------
  # plot : ggplot geom_bar object with error bars.
  
  sig <- data.frame(matrix(nrow = length(levels(plt_df$variable)), ncol=2))
  colnames(sig) <- c('variable', 'adj_p_val')
  for (i in 1:length(levels(plt_df$variable))){
    sig[i,] <- c(levels(plt_df$variable)[i], test_results[test_results$Feature_Name == levels(plt_df$variable)[i],]$adj_p_value)
  }
  sig$adj_p_val <- as.numeric(sig$adj_p_val)
  
  plt_df$high_pt <- plt_df$Mean + plt_df$SE
  if (sum(colnames(plt_df)=='class') > 0){
    max_arr <- acast(plt_df, variable ~ class, value.var = 'high_pt')
  } else {
    max_arr <- acast(plt_df, variable ~ groups, value.var = 'high_pt')
  }
  y <- apply(max_arr, 1, max) + y_shift
  
  sig_marks <- data.frame(matrix(nrow=sum(sig$adj_p_val < alpha), ncol=5))
  colnames(sig_marks) <- c('x', 'xend', 'ylow', 'ylow', 'yhigh')
  j <- 1
  for (i in 1:nrow(sig)){
    if (sig[i,]$adj_p_val < alpha){
      sig_marks[j,] <- c(0.8+(i-1), 1.2+(i-1), y[i], y[i], y[i] + 0.02)
      j <- j + 1
    }
  }
  
  if (nrow(sig_marks) > 0 ){
    for (i in 1:nrow(sig_marks)){
      plot <- plot + sigbra(sig_marks[i,1], sig_marks[i,2], sig_marks[i,3], sig_marks[i,4], sig_marks[i,5], label = '*', text.size = text.size, lab.space = lab.space)
    }
  }
  
  return(plot)
}

plot_class_feature_means <- function(df, df.class, plot_path, ttest_results=NULL, width = NULL, height = NULL, limited = F, prefix = NULL){
  require(reshape2)
  require(ggplot2)
  
  if (limited == T){
    limited_features <- c("total_distance", 
                          "net_distance", 
                          "linearity",
                          "spearmanrsq",
                          "progressivity", 
                          "avg_speed", 
                          "MSD_slope", 
                          "hurst_RS", 
                          "nongauss",
                          "rw_kurtosis01",
                          "rw_kurtosis05",
                          "avg_moving_speed05",
                          "time_moving05",
                          "autocorr_5")
    class_df <- df[limited_features]
    class_df$class <- as.factor(df.class)
  } else {
    class_df <- df
    class_df$class <- as.factor(df.class)
  }
  sem <- function(x){sd(x)/sqrt(length(x))}
  means <- aggregate(. ~ class, data = class_df, FUN = mean)
  stderr <- aggregate(. ~ class, data = class_df, FUN = sem)
  
  means.m <- melt(means, value.name = "Mean")
  stderr.m <- melt(stderr, value.name = "SE")
  plt_df <- merge(means.m, stderr.m)
  write.csv(plt_df, file = paste(plot_path, "class_feature_means.txt"), row.names = F)
  
  if (is.null(width)){
    width = 8
  }
  if (is.null(height)){
    height = 5
  }
  
  feature_theme <- theme(
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1), 
    text= element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )
  
  
  mean_limits <- aes(ymax = plt_df$Mean + plt_df$SE, ymin = plt_df$Mean - plt_df$SE)
  p_mean <- ggplot(plt_df, aes(variable, Mean, fill = class)) + geom_bar(position = "dodge", stat = "identity") + geom_errorbar(mean_limits, position = position_dodge(0.9), width = 0.3) + feature_theme + labs(title = "Feature Means", x = "Feature", y = "Mean (+/- SE)")
  ## add sig bars
  # Find sig levels
  if (!is.null(ttest_results)){
    p_mean <- add_sig_bars(p_mean, plt_df, ttest_results)
  }
  ggsave(p_mean, filename = paste(prefix, 'class_feature_means.png', sep = ''), path = plot_path, width = width, height = height, units = 'in')
  
}

plot_class_feature_medians <- function(df, df.class, plot_path){
  require(reshape2)
  require(ggplot2)
  class_df <- df
  class_df$class <- as.factor(df.class)
  
  sem <- function(x){sd(x)/sqrt(length(x))}
  medians <- aggregate(. ~ class, data = class_df, FUN = median)
  stderr <- aggregate(. ~ class, data = class_df, FUN = sem)
  
  medians.m <- melt(medians, value.name = "Median")
  stderr.m <- melt(stderr, value.name = "SE")
  plt_df <- merge(medians.m, stderr.m)
  write.csv(plt_df, file = paste(plot_path, "class_feature_medians.txt"), row.names = F)
  
  mean_limits <- aes(ymax = plt_df$Median + plt_df$SE, ymin = plt_df$Median - plt_df$SE)
  p_median <- ggplot(plt_df, aes(variable, Median, fill = class)) + geom_bar(position = "dodge", stat = "identity") + geom_errorbar(mean_limits, position = position_dodge(0.9), width = 0.3) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text= element_text(size = 10)) + labs(title = "Feature Means", x = "Feature", y = "Mean (+/- SE)")
  ggsave(p_median, filename = 'feature_medians.png', path = plot_path, width = 8, height = 5, units = 'in')
}

plot_cluster_features <- function(df, df.groups, plot_path, anova_results=NULL, width = 7, height = 4, limited = F, prefix=NULL){
  
  require(reshape2)
  require(ggplot2)
  
  if (limited == T){
    limited_features <- c("total_distance", 
                          "net_distance", 
                          "linearity",
                          "spearmanrsq",
                          "progressivity", 
                          "avg_speed", 
                          "MSD_slope", 
                          "hurst_RS", 
                          "nongauss",
                          "rw_kurtosis01",
                          "rw_kurtosis05",
                          "avg_moving_speed05",
                          "time_moving05",
                          "autocorr_5")
    group_df <- df[limited_features]
    group_df$groups <- as.factor(df.groups)
  } else {
    group_df <- df
    group_df$groups <- as.factor(df.groups)
  }
  
  sem <- function(x){sd(x)/sqrt(length(x))}
  means = aggregate(. ~ groups, data = group_df, FUN = mean)
  medians = aggregate(. ~ groups, data = group_df, FUN = median)
  stderr <- aggregate(. ~ groups, data = group_df, FUN = sem)
  
  # melt data
  means.m <- melt(means, value.name = "Mean")
  medians.m <- melt(medians, value.name = "Median")
  stderr.m <- melt(stderr, value.name = "SE")
  # make data.frame from melted data
  plt_df <- merge(means.m, stderr.m)
  plt_df <- merge(plt_df, medians.m)
  write.csv(plt_df, file = paste(plot_path, prefix, "cluster_feature_means.csv", sep =''), row.names = F)
  
  
  feature_theme <- theme(
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1), 
    text= element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )
  
  mean_limits <- aes(ymax = plt_df$Mean + plt_df$SE, ymin = plt_df$Mean - plt_df$SE)
  p_mean <- ggplot(plt_df, aes(variable, Mean, fill = groups)) + geom_bar(position = "dodge", stat = "identity") + geom_errorbar(mean_limits, position = position_dodge(0.9), width = 0.3) + labs(title = "Feature Means", x = "Feature", y = "Mean (+/- SE)") + feature_theme
  # Add sig bars
  if (!is.null(anova_results)){
    p_mean <- add_sig_bars(p_mean, plt_df, anova_results)
  }
  ggsave(p_mean, filename = paste(prefix, "cluster_feature_means.png", sep =''), path = plot_path, width = width, height = height, units = "in")
  median_limits <- aes(ymax = plt_df$Median + plt_df$SE, ymin = plt_df$Median - plt_df$SE)
  p_median <- ggplot(plt_df, aes(variable, Median, fill = groups)) + geom_bar(position = "dodge", stat = "identity") + geom_errorbar(median_limits, position = position_dodge(0.9), width = 0.3) + labs(title = "Feature Medians", x = "Feature", y = "Median (+/- SE)") + feature_theme
  ggsave(p_median, filename = paste(prefix, "cluster_feature_medians.png", sep =''), path = plot_path, width = width, height = height, units = "in")
  
}

single_class_plots <- function(class_df, plot_path, class_num){
  fgf_sub <- subset(class_df, class_df$class == "FGF2+")
  fgf_sub.noturn <- data.frame(scale(fgf_sub[1:56]))
  nofgf_sub <- subset(class_df, class_df$class == "FGF2-")
  nofgf_sub.noturn <- data.frame(scale(nofgf_sub[1:56]))
  
  fgf_sub.t <- Rtsne(as.matrix(fgf_sub.noturn), perplexity = 50)
  nofgf_sub.t <- Rtsne(as.matrix(nofgf_sub.noturn), perplexity = 50)
  
  fgf_sub.d <- dist(fgf_sub.noturn)
  fgf_sub.fit <- hclust(fgf_sub.d, method = "ward.D")
  fgf_sub.groups <- cutree(fgf_sub.fit, k = class_num)
  
  nofgf_sub.d <- dist(nofgf_sub.noturn)
  nofgf_sub.fit <- hclust(nofgf_sub.d, method = "ward.D")
  nofgf_sub.groups <- cutree(nofgf_sub.fit, k = class_num)
  
  no_axes_bg <- theme(axis.line=element_blank(),
                      axis.text.x=element_blank(),
                      axis.text.y=element_blank(),
                      axis.ticks=element_blank(),
                      axis.title.x=element_blank(),
                      axis.title.y=element_blank(),
                      legend.position="none",
                      panel.background=element_blank(),
                      panel.border=element_blank(),
                      panel.grid.major=element_blank(),
                      panel.grid.minor=element_blank(),
                      plot.background=element_blank())
  
  tsne_groups_fgf <- ggplot(data.frame(fgf_sub.t$Y), aes(x=X1, y=X2, col = as.factor(fgf_sub.groups))) + geom_point(size=1) + labs(title = 'SC FGF+ tSNE') + no_axes_bg
  ggsave(tsne_groups_fgf, path = plot_path, filename = "tsne_groups_fgf_only.png", width = 4, height = 4, units = "in")
  tsne_groups_nofgf <- ggplot(data.frame(nofgf_sub.t$Y), aes(x=X1, y=X2, col = as.factor(nofgf_sub.groups))) + geom_point(size=1) + labs(title = 'SC FGF- tSNE') + no_axes_bg
  ggsave(tsne_groups_nofgf, path = plot_path, filename = "tsne_groups_nofgf_only.png", width = 4, height = 4, units = "in")
  
  plot_cluster_features(fgf_sub.noturn, fgf_sub.groups, file_prefix = "fgfonly_")
  plot_cluster_features(nofgf_sub.noturn, nofgf_sub.groups, file_prefix = "nofgfonly_")
  
  both_subs <- rbind(fgf_sub.noturn, nofgf_sub.noturn)
  both_subs <- data.frame( scale(both_subs) )
  both_subs.groups <- c(fgf_sub.groups, (nofgf_sub.groups + 4))
  
  plot_cluster_features(both_subs, both_subs.groups, file_prefix = "bothsub_")
}

## Statistics

fgf_ttest_feature <- function(class_df, feature_name, plot_path, prefix=''){
  fgf2_sub <- subset(class_df, class_df$class == "FGF2+")
  nofgf2_sub <- subset(class_df, class_df$class == "FGF2-")
  
  t_result <- t.test(fgf2_sub[,feature_name], nofgf2_sub[,feature_name], alternative = "two.sided")
  capture.output(t_result, file = paste(plot_path, prefix, "ttest_", feature_name, ".txt", sep=''))
  return(t_result$p.value)
}

mef_ttest_feature <- function(class_df, feature_name, plot_path, prefix=''){
  myc_sub <- subset(class_df, class_df$class == "MycRas")
  wt_sub <- subset(class_df, class_df$class == "WT")
  
  t_result <- t.test(myc_sub[,feature_name], wt_sub[,feature_name], alternative = "two.sided")
  capture.output(t_result, file = paste(plot_path, prefix, "ttest_", feature_name, ".txt", sep=''))
  return(t_result$p.value)
}

anova_feature <- function(df, df.groups, feature_name, plot_path, prefix=''){
  groups <- as.factor(df.groups)
  av <- aov(df[,feature_name] ~ groups, data = df)
  capture.output(summary(av), file = paste(plot_path, prefix, "anova_", feature_name, ".txt", sep =''))
  p <- summary(av)[[1]][['Pr(>F)']][1]
  return(p)
}

fgf_contingency_table_preftest <- function(df.cluster_membership, df.groups, plot_path, experiment = NULL){
  ct <- matrix(0, ncol = max(df.groups), nrow = 2)
  ct[1,] <- df.cluster_membership$fgf
  ct[2,] <- df.cluster_membership$nofgf
  
  # Fisher test
  fisher <- fisher.test(ct)
  capture.output(fisher, file = paste(plot_path, experiment, '_fisher.txt', sep = ''))
  # Chi^2 test
  chi2 <- chisq.test(ct)
  capture.output(fisher, file = paste(plot_path, experiment, '_chi2.txt', sep = ''))
  
}

mef_contingency_table_preftest <- function(df.cluster_membership, df.groups, plot_path, experiment = NULL){
  ct <- matrix(0, ncol = max(df.groups), nrow = 2)
  ct[1,] <- df.cluster_membership$mycras
  ct[2,] <- df.cluster_membership$wt
  
  # Fisher test
  fisher <- fisher.test(ct)
  capture.output(fisher, file = paste(plot_path, experiment, '_fisher.txt', sep = ''))
  # Chi^2 test
  chi2 <- chisq.test(ct)
  capture.output(chi2, file = paste(plot_path, experiment, '_chi2.txt', sep = ''))
  
}

density_contingency_table_preftest <- function(groups_mat, plot_path, prefix=''){
  chi2 <- chisq.test(t(groups_mat))
  capture.output(chi2, file = paste(plot_path, prefix, '_chi2.txt', sep = ''))
}

density_regression_tests <- function(df, density){
  df$density <- density
  
  features <- colnames(df)
  features <- features[features != 'density']
  
  results <- data.frame(matrix(nrow=length(features), ncol=3, 0))
  colnames(results) <- c('Fstat', 'Rsq', 'p')
  
  for (i in 1:length(features)){
    f <- features[i]
    f.lm <- lm(df[,f] ~ density, df)
    f.s <- summary(f.lm)
    rsq <- f.s$r.squared
    fstat <- f.s$fstatistic[1]
    p <- pf(f.s$fstatistic[1], f.s$fstatistic[2], f.s$fstatistic[3], lower.tail = F)
    results[i,] <- c(fstat, rsq, p)
  }
  results$feature <- features
  return(results)
}

#####

real_units_speed <- function(df.raw, df.class, time_interval=6.5, unit_multiplier=0.325){
  # Find the mean displacement speed of each cell category in df.class
  # in real units, as specified by a multiplier in unit_multiplier
  # which converts pixels to real distance
  # i.e. 6.5 um / px on the Hamamatsu camera
  # Parameters
  # ----------
  # df.raw : data frame with raw (unscaled) feature values
  # df.class : vector of class labels
  # time_interval : float, time in minutes per frame.
  # unit_multiplier : float, um distance per pixel (make sure to include effect of Mag!).
  # Returns
  # -------
  # real_speeds : data frame, k x 3, class name, <speed in um/min>, speed SEM 
  #
  
  real_speeds = data.frame(matrix(nrow = length(unique(df.class)), ncol = 3))
  for ( i in 1:length(unique(df.class)) ){
    real_speeds[i,] = c(unique(df.class)[i], 
                        mean(df.raw[df.class==unique(df.class)[i],]$avg_speed)/time_interval * unit_multiplier,
                        (sd(df.raw[df.class==unique(df.class)[i],]$avg_speed)/time_interval * unit_multiplier)/sqrt(sum(df.class==unique(df.class)[i])))
  }
  colnames(real_speeds) <- c('Class', 'SpeedMean', 'SpeedSE')
  return(real_speeds)
}

## Cluster Validation

clustval <- function(df, class_num, plot_path, prefix = NULL, method = 'ward.D2'){
  df.groups <- cutree(hclust(dist(df), method = "ward.D"), k = class_num)  
  
  require(NbClust)
  df.nbclust <- NbClust(df, min.nc = 2, max.nc = (class_num + 3), method = method)
  capture.output(df.nbclust$All.index, file = paste(plot_path, prefix, "nbclust_indices.txt", sep = ''))
  capture.output(df.nbclust$All.CriticalValues, file = paste(plot_path, prefix, "nbclust_critvals.txt", sep = ''))
  capture.output(df.nbclust$Best.nc, file = paste(plot_path, prefix, "nbclust_bestnc.txt", sep = ''))
  
  # library(clv)
  # df.scattd <- cls.scatt.data(df, df.groups, dist="euclidean")
  # df.conn <- connectivity(df, df.groups, neighbour.num = 10)
  # df.dunn <- clv.Dunn(df.scattd, "complete", intercls = "complete")
  
  #library(clusteval)
  #hclust_fun <- function(x, num_clusters){ cutree(hclust(dist(x), method = "ward.D"), num_clusters)}
  #df.clustomit <- clustomit(df, class_num, hclust_fun)
  #capture.output(df.clustomit$boot_similarity, file = paste(plot_path, prefix, "clustomit_bootsimilarity.txt", sep = ''))
  #capture.output(df.clustomit$boot_aggregate, file = paste(plot_path, prefix, "clustomit_bootaggregate.txt", sep = ''))
  
}

### LESS COMMON


## Cluster centroid distances

centroid_distances <- function(df, df.groups){
  # Takes data.frame and factor of cluster assignments
  # Returns distance matrix b/w centroids of each cluster
  # d = [d(i1,i2), ]
  #     [d(i1,i3), d(i2,i3)], &c.
  numc <- max(df.groups)
  centroids <- matrix(0, numc, ncol(df))
  for (i in 1:numc){
    clust <- df[which(df.groups == i),]
    centroid <- colMeans(clust, na.rm = T)
    centroids[i,] = centroid
  }
  d <- dist(centroids, diag = T, upper = T)
  
  return(d)
}

diff_centroid_distance <- function(class_df, df.groups){
  fgf2_sub.groups <- df.groups[which(class_df$class == "FGF2+")]
  nofgf2_sub.groups <- df.groups[which(class_df$class == "FGF2-")]
  fgf2_sub <- subset(class_df, class_df$class == "FGF2+")
  nofgf2_sub <- subset(class_df, class_df$class == "FGF2-")
  fgf2_sub$class <- NULL
  nofgf2_sub$class <- NULL
  
  fgf2_sub.d <- centroid_distances(fgf2_sub, fgf2_sub.groups)
  nofgf2_sub.d <- centroid_distances(nofgf2_sub, nofgf2_sub.groups)
  
  diff_d <- fgf2_sub.d - nofgf2_sub.d
  
  return(diff_d)
}

plot_diff_cent_dist <- function(centroid_dist, plot_path, prefix = NULL){
  require(reshape2)
  require(ggplot2)
  d.m <- melt(as.matrix(centroid_dist))
  diff_theme <- theme(#legend.position = "none",
    #axis.text.x=element_blank(),
    #axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    legend.position="right",
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank())
  p <- ggplot(d.m, aes(x = Var2, y = Var1, fill = value))
  p <- p + geom_tile() + labs(x = "Cluster", y = "Cluster") + diff_theme + coord_fixed() + scale_fill_distiller(palette = "RdBu", limits=c(-1.5,1.5))
  ggsave(p, filename= paste(prefix, "centroid_distance_diff.png", sep = ''), path = plot_path, width = 4, height = 4)
}


## ICA

ica_plots <- function(df, plot_path, class_num){
  require(ica)
  require(ggplot2)
  df.ica <- icafast(df, n = class_num, center = TRUE, maxit = 1000, alpha = 1)
  df.ica.d <- dist(df.ica$S, method = "euclidean")
  df.ica.fit <- hclust(df.ica.d, method = 'ward.D')
  df.ica.groups <- as.factor(cutree(df.ica.fit, k = class_num))
  ic_labs <- c("IC1", "IC2", "IC3", "IC4", "IC5", "IC6")
  df.ics <- data.frame(df.ica$S)
  colnames(df.ics) <- ic_labs[1:ncol(df.ics)]
  ica_theme <- theme(legend.position = "none",
                     panel.background=element_blank(),
                     panel.border=element_blank(),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),
                     plot.background=element_blank(),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks = element_blank())
  ica_wardD <- ggplot(df.ics, aes(x = IC1, y = IC2, col = df.ica.groups)) + geom_point() + ica_theme + labs(title= "ICA" )
  ggsave(ica_wardD, filename = "ica_wardD.png", path = plot_path, width = 4, height = 4, units = "in")
}

sem <- function(x){sd(x, na.rm = T)/sqrt(length(x) - sum(is.nan(x)))}

