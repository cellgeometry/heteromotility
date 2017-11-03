# Myoblast Multiexperiment Analysis

source('hm_common_fnx.R')


process_myoblasts <- function(myo_exp1, myo_exp2, plot_path, class_num, method='ward.D2', nb_clust=F, tsne_run=F){
  
  plot_path <- "data/myoblast/"
  
  fgf_df1 <- read.csv(paste('data/', myo_exp1, '/', 'fgf2_exp_motility_statistics.csv', sep = ''))
  nofgf_df1 <- read.csv(paste('data/', myo_exp1, '/', 'nofgf2_exp_motility_statistics.csv', sep = ''))
  fgf_df2 <- read.csv(paste('data/', myo_exp2, '/', 'fgf2_exp_motility_statistics.csv', sep = ''))
  nofgf_df2 <- read.csv(paste('data/', myo_exp2, '/', 'nofgf2_exp_motility_statistics.csv', sep = ''))
  
  plate1 = rbind(fgf_df1, nofgf_df1)
  plate1.class <- c(rep("FGF2+", nrow(fgf_df1)), rep("FGF2-", nrow(nofgf_df1)))
  plate1.class_num <- c(rep(1, nrow(fgf_df1)), rep(2, nrow(nofgf_df1)))
  plate2 = rbind(fgf_df2, nofgf_df2)
  plate2.class <- c(rep("FGF2+", nrow(fgf_df2)), rep("FGF2-", nrow(nofgf_df2)))
  plate2.class_num <- c(rep(1, nrow(fgf_df2)), rep(2, nrow(nofgf_df2)))
  
  r1 = c(6,15,19,43,48,49,51,53,66,70,115,135,150,154,167,202,204,205,208)
  r2 = c(3, 18, 19, 92, 97, 104, 112)
  
  plate1 = plate1[-r1,]
  plate1.class = plate1.class[-r1]
  plate1.class_num = plate1.class_num[-r1]
  plate2 = plate2[-r2,]
  plate2.class = plate2.class[-r2]
  plate2.class_num = plate2.class_num[-r2]
  
  df <- rbind(plate1, plate2)
  df.plate <- c( rep(1, nrow(plate1)), rep(2, nrow(plate2)) )
  df.class <- c(plate1.class, plate2.class)
  df.class_num <- c(plate1.class_num, plate2.class_num)
  df <- df[ -c(1,2) ]
  df.raw <- df
  df <- data.frame( scale (df) )
  df.plate <- df.plate[!duplicated(df)]
  df.raw <- df.raw[!duplicated(df),]
  df.class <- df.class[!duplicated(df)]
  df.class_num <- df.class_num[!duplicated(df)]
  df <- df[!duplicated(df),]
  df.noturn <- df[,1:55]
  class_df.noturn <- df.noturn
  class_df.noturn$class <- df.class
  class_df.raw <- df.raw
  class_df.raw$class <- df.class
  
  pca_loading_variation(df.noturn, plot_path, prefix='Myoblast_')
  
  df.comp <- pca_plots(df.noturn, df.class, df.plate, plot_path, class_num)
  #df.g <- heatmap3_plots(df.noturn, df.class=df.class, plot_path, class_num, method=method)
  df.groups <- cutree(hclust(dist(df.comp[,1:30]), method=method), class_num)
  if (tsne_run){tsne_plots(df.comp, df.groups, df.class, df.plate, plot_path, width=2.5, height=2.5)}
  #som_plots(df.comp, grid_size = 20, plot_path = plot_path, class_num)
  
  # plot clustograms with different linkage metrics
  for (m in c('ward.D2', 'average', 'complete', 'median')){
    heatmap2_plots(df.noturn, df.class = df.class, plot_path, class_num, method=m, prefix=paste('clustogram_', m, '_', sep=''))
  }
  
  ttest_results = data.frame(matrix(ncol=3, nrow=ncol(df.raw)))
  for (i in 1:ncol(df.raw)){
    feature_name <- colnames(df.raw)[i]
    p <- fgf_ttest_feature(class_df.raw, feature_name, paste(plot_path, "ttests/", sep =''))
    ttest_results[i,] <- c(feature_name, p, 0)
  }
  colnames(ttest_results) <- c('Feature_Name', 'p_value', 'adj_p_value')
  ttest_results$adj_p_value <- p.adjust(ttest_results$p_value, method = 'holm')
  ttest_results <- ttest_results[order(ttest_results$adj_p_value),] # order by pvalue, ascending
  write.csv(ttest_results, file = paste(plot_path, 'ttests/', '0_ttest_summary.csv', sep=''), row.names = F)
  
  anova_results = data.frame(matrix(ncol=3, nrow=ncol(df.raw)))
  for (i in 1:ncol(df.raw)){
    feature_name <- colnames(df.raw)[i]
    p <- anova_feature(df.raw, df.groups, feature_name, paste(plot_path, "anovas/", sep = ''))
    anova_results[i,] <- c(feature_name, p, 0)
  }
  colnames(anova_results) <- c('Feature_Name', 'p_value', 'adj_p_value')
  anova_results$adj_p_value <- p.adjust(anova_results$p_value, method = 'holm')
  anova_results <- anova_results[order(anova_results$adj_p_value),] # order by pvalue, ascending
  write.csv(anova_results, file = paste(plot_path, 'anovas/', '0_anova_summary.csv', sep=''), row.names = F)
  
  plot_class_feature_means(df.noturn, df.class, plot_path, ttest_results=ttest_results, width = 7, height = 4)
  plot_class_feature_means(df.noturn, df.class, plot_path, ttest_results=ttest_results, width = 4, height = 3, limited = T, prefix = "limited_")
  plot_class_feature_medians(df.noturn, df.class, plot_path)
  plot_cluster_features(df, df.groups, plot_path, anova_results=anova_results, width = 7, height = 4)
  plot_cluster_features(df, df.groups, plot_path, anova_results=anova_results, width = 4, height = 3, limited = T, prefix = "limited_")
  
  df.cluster_membership <- fgf_cluster_props(class_df.noturn, df.groups)
  fgf_plot_cluster_memberships(df.cluster_membership, plot_path)
  fgf_contingency_table_preftest(df.cluster_membership, df.groups, plot_path, experiment = "Myoblast_ClustPref")
  #centroid_distance_shift <- diff_centroid_distance(class_df, df.groups)
  #plot_diff_cent_dist(centroid_distance_shift, plot_path)
  
  
  classif_fgf_df <- df.noturn
  classif_fgf_df$class <- df.class_num
  classif_groups_df <- df.noturn
  classif_groups_df$group <- df.groups
  write.csv(classif_fgf_df, file = paste(plot_path, "myoblast_fgf_classif.csv", sep =''), row.names = F)
  write.csv(classif_groups_df, file = paste(plot_path, "myoblast_cluster_classif.csv", sep =''), row.names = F)
  
  # Save speed in real units
  real_speeds = real_units_speed(df.raw, df.class, time_interval = 6.5, unit_multiplier = 0.325)
  write.csv(real_speeds, file = paste(plot_path, 'Myoblast_speeds_real_units.csv', sep=''), row.names = F)
  
  mv <- manova(as.matrix(df.comp) ~ as.factor(df.groups))
  capture.output(summary(mv), file = paste(plot_path, "manova_summary.txt", sep = ''))
  
  if (nb_clust == TRUE){
    clustval(df.comp, class_num = class_num, plot_path = plot_path, prefix = "Myoblast_ClustVal_", method=method)
  }
}

myo_exp1 <- "myoblast/20160623"
myo_exp2 <- "myoblast/20160720"
method = 'ward.D2'
class_num = 2
nb_clust = T
tsne_run = T
process_myoblasts(myo_exp1, myo_exp2, plot_path, class_num=class_num, method=method, nb_clust=nb_clust, tsne_run = tsne_run)

plot_path = 'data/myoblast/resampling'
for (class_num in 3:4){
  for (i in 1:10){
    idx <- sample(1:nrow(df.comp), size = as.integer(nrow(df.comp)*0.8))
    df <- df.comp[idx,]
    
    df.groups <- cutree(hclust(dist(df[,1:30], method='euclidean'), method=method), class_num)
    p <- ggplot(df, aes(x=PC1, y=PC2, color=as.factor(df.groups))) + geom_point() + guides(fill=guide_legend(title="Cluster"))
    ggsave(p, path = plot_path, filename=paste('pca_groups_k', class_num, '_iter', i, '.png', sep=''), width=4, height=3, units='in')
  }
}
