# Plate to plate variation 

source('hm_common_fnx.R')

compare_plates <- function(experiment1, experiment2, experiment3, class_num, method='ward.D', nb_clust=F, tsne_run=F){
  
  plot_path = "data/musc/"
  
  fgf_df1 = read.csv( paste('data/', experiment1, '/', 'fgf2_exp_motility_statistics.csv', sep = '') )
  nofgf_df1 = read.csv( paste('data/', experiment1, '/', 'nofgf2_exp_motility_statistics.csv', sep = '') )
  
  fgf_df2 = read.csv( paste('data/', experiment2, '/', 'fgf2_exp_motility_statistics.csv', sep = '') )
  nofgf_df2 = read.csv( paste('data/', experiment2, '/', 'nofgf2_exp_motility_statistics.csv', sep = '') )
  
  fgf_df3 = read.csv( paste('data/', experiment3, '/', 'fgf2_exp_motility_statistics.csv', sep = '') )
  nofgf_df3 = read.csv( paste('data/', experiment3, '/', 'nofgf2_exp_motility_statistics.csv', sep = '') )
  
  plate1 = rbind(fgf_df1, nofgf_df1)
  plate1.class <- c(rep("FGF2+", nrow(fgf_df1)), rep("FGF2-", nrow(nofgf_df1)))
  plate1.class_num <- c(rep(1, nrow(fgf_df1)), rep(2, nrow(nofgf_df1)))
  plate2 = rbind(fgf_df2, nofgf_df2)
  plate2.class <- c(rep("FGF2+", nrow(fgf_df2)), rep("FGF2-", nrow(nofgf_df2)))
  plate2.class_num <- c(rep(1, nrow(fgf_df2)), rep(2, nrow(nofgf_df2)))
  plate3 = rbind(fgf_df3, nofgf_df3)
  plate3.class <- c(rep("FGF2+", nrow(fgf_df3)), rep("FGF2-", nrow(nofgf_df3)))
  plate3.class_num <- c(rep(1, nrow(fgf_df3)), rep(2, nrow(nofgf_df3)))
  
  r <- c(49, 84, 624, 1238, 1386, 1472, 1599, 1738, 1776, 1796, 1907, 1954, 1963,  2091, 2468, 2543, 2989, 3624, 3873, 4006, 4145, 4192)
  
  all_plates <- rbind(plate1, plate2, plate3)
  all_plates.class <- c(plate1.class, plate2.class, plate3.class)
  all_plates.class_num <- c(plate1.class_num, plate2.class_num, plate3.class_num)
  
  all_plates <- all_plates[-r,]
  all_plates.class <- all_plates.class[-r]
  all_plates.class_num <- all_plates.class_num[-r]
  
  all_plates <- all_plates[ -c(1,2) ]
  all_plates.plate <- c( rep(1, nrow(plate1)), rep(2, nrow(plate2)), rep(3, nrow(plate3)))
  all_plates.plate <- all_plates.plate[-r]
  all_plates.raw <- all_plates
  all_plates <- data.frame( scale(all_plates) )
  all_plates.raw <- all_plates.raw[!duplicated(all_plates),]
  all_plates.plate <- all_plates.plate[!duplicated(all_plates)]
  all_plates.class <- all_plates.class[!duplicated(all_plates)]
  all_plates.class_num <- all_plates.class_num[!duplicated(all_plates)]
  all_plates <- all_plates[!duplicated(all_plates),]
  all_plates.noturn <- all_plates[,1:55]
  class_all_plates <- all_plates
  class_all_plates$class <- all_plates.class
  class_all_plates.noturn <- all_plates.noturn
  class_all_plates.noturn$class <- all_plates.class
  class_all_plates.raw <- all_plates.raw
  class_all_plates.raw$class <- all_plates.class
  
  
  all_plates.comp <- pca_plots(all_plates.noturn, all_plates.class, all_plates.plate, plot_path, class_num)
  all_plates.groups <- cutree(hclust(dist(all_plates.comp[,1:30], method='euclidean'), method=method), class_num)

  #all_plates.g <- heatmap3_plots(all_plates.noturn, df.class=all_plates.class, plot_path, class_num, method = method)
  if (tsne_run){tsne_plots(all_plates.comp, all_plates.groups, all_plates.class, all_plates.plate, plot_path, perplexity=70)}
  #tsne_single_class_plots(class_all_plates, plot_path, class_num)
  #som_plots(all_plates.comp, grid_size = 20, plot_path = plot_path, class_num)
  
  pca_loading_variation(all_plates.noturn, plot_path, prefix='MuSC_')
  
  ttest_results = data.frame(matrix(ncol=3, nrow=ncol(all_plates.raw)))
  for (i in 1:ncol(all_plates.raw)){
    feature_name <- colnames(all_plates.raw)[i]
    p <- fgf_ttest_feature(class_all_plates.raw, feature_name, paste(plot_path, "ttests/", sep =''))
    ttest_results[i,] <- c(feature_name, p, 0)
  }
  colnames(ttest_results) <- c('Feature_Name', 'p_value', 'adj_p_value')
  ttest_results$adj_p_value <- p.adjust(ttest_results$p_value, method = 'holm')
  ttest_results <- ttest_results[order(ttest_results$adj_p_value),] # order by pvalue, ascending
  write.csv(ttest_results, file = paste(plot_path, 'ttests/', '0_ttest_summary.csv', sep=''), row.names = F)
  
  anova_results = data.frame(matrix(ncol=3, nrow=ncol(all_plates.raw)))
  for (i in 1:ncol(all_plates.raw)){
    feature_name <- colnames(all_plates.raw)[i]
    p <- anova_feature(all_plates.raw, all_plates.groups, feature_name, paste(plot_path, "anovas/", sep = ''))
    anova_results[i,] <- c(feature_name, p, 0)
  }
  colnames(anova_results) <- c('Feature_Name', 'p_value', 'adj_p_value')
  anova_results$adj_p_value <- p.adjust(anova_results$p_value, method = 'holm')
  anova_results <- anova_results[order(anova_results$adj_p_value),] # order by pvalue, ascending
  write.csv(anova_results, file = paste(plot_path, 'anovas/', '0_anova_summary.csv', sep=''), row.names = F)
  
  plot_class_feature_means(all_plates.noturn, all_plates.class, plot_path, ttest_results = ttest_results, width = 7, height = 4)
  plot_class_feature_means(all_plates.noturn, all_plates.class, plot_path, ttest_results = ttest_results, width = 4, height = 3, limited = T, prefix = "limited_")
  plot_class_feature_medians(all_plates.noturn, all_plates.class, plot_path)
  plot_cluster_features(all_plates, all_plates.groups, plot_path, anova_results = anova_results, width = 7, height = 4)
  plot_cluster_features(all_plates, all_plates.groups, plot_path, anova_results = anova_results, width = 4, height = 3, limited = T, prefix = "limited_")
  
  df.cluster_membership <- fgf_cluster_props(class_all_plates.noturn, all_plates.groups)
  fgf_plot_cluster_memberships(df.cluster_membership, plot_path)
  fgf_contingency_table_preftest(df.cluster_membership, all_plates.groups, plot_path, experiment = "MuSC_ClustPref")
  centroid_distance_shift <- diff_centroid_distance(class_all_plates, all_plates.groups)
  plot_diff_cent_dist(centroid_distance_shift, plot_path)
  
  
  classif_fgf_df <- all_plates.noturn
  classif_fgf_df$class <- all_plates.class_num
  classif_groups_df <- all_plates.noturn
  classif_groups_df$group <- all_plates.groups
  write.csv(classif_fgf_df, file = paste(plot_path, "sc_fgf_classif.csv", sep =''), row.names = F)
  write.csv(classif_groups_df, file = paste(plot_path, "sc_cluster_classif.csv", sep =''), row.names = F)
  
  # Save speeds in real units
  real_speeds = real_units_speed(all_plates.raw, all_plates.class, time_interval = 6.5, unit_multiplier = 0.325)
  write.csv(real_speeds, file = paste(plot_path, 'MuSC_speeds_real_units.csv', sep=''), row.names = F)
  real_speeds_clusters = real_units_speed(all_plates.raw, all_plates.groups, time_interval = 6.5, unit_multiplier = 0.325)
  write.csv(real_speeds_clusters, file = paste(plot_path, 'MuSC_speeds_real_units_clusters.csv', sep=''), row.names = F)
  
  
  mv <- manova(as.matrix(all_plates.comp) ~ as.factor(all_plates.groups))
  capture.output(summary(mv), file = paste(plot_path, "manova_summary.txt", sep = ''))
  
  if (nb_clust == TRUE){
    clustval(all_plates.comp, class_num, plot_path = plot_path, prefix = paste("MuSC_ClustVal_", method, '_', sep = ''), method = method)
    #clustval(all_plates.comp, class_num, plot_path = plot_path, prefix = paste("MuSC_ClustVal_", 'KMeans', '_', sep = ''), method = 'kmeans')
  }
}

## Execute

experiment1 = 'musc/20160626'
experiment2 = 'musc/20160701_0'
experiment3 = 'musc/20160701_1'
class_num = 4
method = 'ward.D'
nb_clust = F
tsne_run = F
compare_plates(experiment1, experiment2, experiment3, class_num, method = method, nb_clust=nb_clust)


# resampling

plot_path = 'data/musc/resampling/'


for (class_num in 3:4){
  for (i in 1:10){
    idx <- sample(1:nrow(all_plates.comp), size = as.integer(nrow(all_plates.comp)*0.8))
    df <- all_plates.comp[idx,]
    
    df.groups <- cutree(hclust(dist(df[,1:30], method='euclidean'), method=method), class_num)
    p <- ggplot(df, aes(x=PC1, y=PC2, color=as.factor(df.groups))) + geom_point() + guides(fill=guide_legend(title="Cluster"))
    ggsave(p, path = plot_path, filename=paste('pca_groups_k', class_num, '_iter', i, '.png', sep=''), width=4, height=3, units='in')
  }
}
