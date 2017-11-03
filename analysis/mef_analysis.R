## MycRas v WT Plots

source('hm_common_fnx.R')

mycras_wt_plots <- function(myc_exps, wt_exps, class_num, method='ward.D2', nb_clust=F, tsne_run=F, plot_path='data/mef/', prefix=''){
  # Runs analysis and generates plots for MycRas and WT MEF data. 
  #
  # Parameters
  # ----------
  # XYZ_exp : string used to name first exp.
  # class_num : integer. number of classes to specify for functions with a k parameter.
  # method : string. method argument for hierarchical clustering. see hclust docs for options.
  # nb_clust : boolean. Run NbClust.
  # tsne_run : boolean. Run tsne and produce plots.
  #
  # Returns
  # -------
  # None.
  #

  # Load experiments from lists of experiment string locations
  myc_exp_sizes = c()
  myc_df <- data.frame()
  for (mexp in myc_exps){
    m_df = read.csv(paste('data/', mexp, '/', 'exp_motility_statistics.csv', sep =''))
    myc_exp_sizes = c(myc_exp_sizes, nrow(m_df))
    myc_df = rbind(myc_df, m_df)
  }
  
  wt_exp_sizes = c()
  wt_df <- data.frame()
  for (wtexp in wt_exps){
    wexp_df = read.csv(paste('data/', wtexp, '/', 'exp_motility_statistics.csv', sep =''))
    wt_exp_sizes = c(wt_exp_sizes, nrow(wexp_df))
    wt_df = rbind(wt_df, wexp_df)
  }
  all_exp_sizes <- c(myc_exp_sizes, wt_exp_sizes)
  df <- rbind(myc_df, wt_df)
  
  # myc_df <- read.csv(paste('data/', myc_exp, '/', 'exp_motility_statistics.csv', sep =''))
  # myc2_df <- read.csv(paste('data/', myc2_exp, '/', 'exp_motility_statistics.csv', sep =''))
  # wt_df <- read.csv( paste('data/', wt_exp, '/', 'exp_motility_statistics.csv', sep ='') )
  # wt2_df <- read.csv( paste('data/', wt2_exp, '/', 'exp_motility_statistics.csv', sep = ''))
  # wt3_df <- read.csv( paste('data/', wt3_exp, '/', 'exp_motility_statistics.csv', sep = ''))
  # wt4_df <- read.csv( paste('data/', wt4_exp, '/', 'exp_motility_statistics.csv', sep = ''))
  # df <- rbind(myc_df, myc2_df, wt_df, wt2_df, wt3_df, wt4_df)
  df <- df[ -c(1,2) ]
  df.raw <- df
  df <- data.frame( scale(df) )
  df.noturn <- df[,1:55]
  # df.class <- c( rep("MycRas", (nrow(myc_df) + nrow(myc2_df) )), rep("WT", (nrow(wt_df) + nrow(wt2_df) + nrow(wt3_df) + nrow(wt4_df)) ) )
  # df.class_num <- c( rep(1, (nrow(myc_df) + nrow(myc2_df) )), rep(2, (nrow(wt_df) + nrow(wt2_df) + nrow(wt3_df) + nrow(wt4_df)) ) )
  # df.plate <- c(rep(1, nrow(myc_df)), rep(2, nrow(myc2_df)), rep(3, nrow(wt_df)), rep(4, nrow(wt2_df)), rep(5, nrow(wt3_df)), rep(6, nrow(wt4_df)))
  df.class <- c(rep("MycRas", nrow(myc_df)), rep('WT', nrow(wt_df)))
  df.class_num <- c(rep(1, nrow(myc_df)), rep(2, nrow(wt_df)))
  df.plate <- c()
  for (i in 1:length(all_exp_sizes)){
    df.plate <- c(df.plate, rep(i, all_exp_sizes[i]))
  }
  class_df <- df
  class_df$class <- as.factor(df.class)
  class_df.raw <- df.raw
  class_df.raw$class <- as.factor(df.class)
  
  pca_loading_variation(df.noturn, plot_path, prefix='MEF_')
  
  df.comp <- pca_plots(df.noturn, df.class, df.plate, plot_path, class_num, prefix=prefix)
  df.g <- heatmap3_plots(df.noturn, df.class, plot_path, class_num, method = method, prefix=prefix)
  df.groups <- cutree(hclust(dist(df.comp[,1:30]), method = method), class_num)
  if (tsne_run){tsne_plots(df.noturn, df.groups, df.class, df.plate, plot_path, method = method, width = 2, height = 2, experiment = prefix)}
  #som_plots(df.noturn, 20, class_num, plot_path, method = method)
  
  ttest_results = data.frame(matrix(ncol=3, nrow=ncol(df.raw)))
  for (i in 1:ncol(df.raw)){
    feature_name <- colnames(df.raw)[i]
    p <- mef_ttest_feature(class_df.raw, feature_name, paste(plot_path, "ttests/", sep =''), prefix=prefix)
    ttest_results[i,] <- c(feature_name, p, 0)
  }
  colnames(ttest_results) <- c('Feature_Name', 'p_value', 'adj_p_value')
  ttest_results$adj_p_value <- p.adjust(ttest_results$p_value, method = 'holm')
  ttest_results <- ttest_results[order(ttest_results$adj_p_value),] # order by pvalue, ascending
  write.csv(ttest_results, file = paste(plot_path, 'ttests/', '0_ttest_summary.csv', sep=''), row.names = F)
  
  anova_results = data.frame(matrix(ncol=3, nrow=ncol(df.raw)))
  for (i in 1:ncol(df.raw)){
    feature_name <- colnames(df.raw)[i]
    p <- anova_feature(df.raw, df.groups, feature_name, paste(plot_path, "anovas/", sep = ''), prefix=prefix)
    anova_results[i,] <- c(feature_name, p, 0)
  }
  colnames(anova_results) <- c('Feature_Name', 'p_value', 'adj_p_value')
  anova_results$adj_p_value <- p.adjust(anova_results$p_value, method = 'holm')
  anova_results <- anova_results[order(anova_results$adj_p_value),] # order by pvalue, ascending
  write.csv(anova_results, file = paste(plot_path, 'anovas/', '0_anova_summary.csv', sep=''), row.names = F)
  
  plot_class_feature_means(df.noturn, df.class, plot_path, ttest_results=ttest_results, width = 7, height = 4, prefix=prefix)
  plot_class_feature_means(df.noturn, df.class, plot_path, ttest_results=ttest_results, width = 4, height = 3, limited = T, prefix = paste(prefix, "limited_", sep=''))
  plot_class_feature_medians(df, df.class, plot_path)
  df.cluster_membership <- mef_cluster_props(class_df, df.groups)
  mef_plot_cluster_memberships(df.cluster_membership, plot_path, prefix=prefix)
  mef_contingency_table_preftest(df.cluster_membership, df.groups, plot_path, experiment = paste(prefix, "MEF_ClustPref", sep=''))
  plot_cluster_features(df.noturn, df.groups, plot_path, anova_results = anova_results, width = 7, height = 4, prefix=prefix)
  plot_cluster_features(df.noturn, df.groups, plot_path, anova_results = anova_results, width = 4, height = 3, limited = T, prefix = paste(prefix, "limited_", sep=''))
  
  capture.output(df.cluster_membership, file = paste(plot_path, prefix, "cluster_membership.txt", sep = ''))
  
  classif_df <- df.noturn
  classif_df$class <- df.class_num
  write.csv(classif_df, file = paste(plot_path, prefix, "mycwt_classification_df.csv", sep =''), row.names = F)
  classif_df$plate <- df.plate
  write.csv(classif_df, file = paste(plot_path, prefix, "mycwt_classification_plate_df.csv", sep =''), row.names = F)
  
  
  # Save speed in real units
  real_speeds = real_units_speed(df.raw, df.class, time_interval = 6.5, unit_multiplier = 0.325)
  write.csv(real_speeds, file = paste(plot_path, prefix, 'MEF_speeds_real_units.csv', sep=''), row.names = F)
  
  cluster_classif_df <- df.noturn
  cluster_classif_df$class <- df.groups
  write.csv(cluster_classif_df, file = paste(plot_path, prefix, "mef_cluster_classification_df.csv", sep =''), row.names = F)
  
  if (nb_clust == TRUE){
    clustval(df.comp, class_num = class_num, plot_path = plot_path, prefix = paste(prefix, "MEF_ClustVal_", sep=''), method=method)
  }
}

myc_exp = "mef/mycras/20160917"
myc2_exp = "mef/mycras/20160918"
wt_exp = "mef/wt/20160711"
wt2_exp = "mef/wt/20160925_0"
wt3_exp = "mef/wt/20160925_1"
wt4_exp = "mef/wt/20160927"
class_num = 3
method='ward.D2'
nb_clust = F
tsne_run = F

myc_exps <- c(myc_exp, myc2_exp)
wt_exps <- c(wt_exp, wt2_exp, wt3_exp, wt4_exp)

# Run all exps together
mycras_wt_plots(myc_exps, wt_exps, class_num=class_num, method=method, nb_clust=nb_clust, tsne_run=tsne_run, prefix='test_')


## Run a round robin for mycras v wt plots
wt_comb <- combn(wt_exps, 3)

wt_grps <- c()
grps_hist <- c()
for (i in 1:ncol(wt_comb)){
 wt_grp <- wt_comb[,i]
 myc_grp <- myc_exps[sample(1:2, 1)]
 tsne_run = T
 grps_hist <- rbind(grps_hist, c(myc_grp, wt_grp))
 mycras_wt_plots(myc_grp, wt_grp, class_num=class_num, method=method, nb_clust=nb_clust, tsne_run=tsne_run, prefix=paste('rrobin_', sprintf('%02d', i), '_', sep=''), plot_path='data/mef/round_robin/')
 
}
write.csv(grps_hist, file = 'data/mef/round_robin/00_groups_index.csv')

## Perform clustering with limited feature subsets 


mycras_wt_feature_reduction <- function(myc_exps, wt_exps, class_num, method='ward.D2', nb_clust=F, tsne_run=F, plot_path='data/mef/feat_red/', prefix='', min_feat_set=F){
  # Runs analysis and generates plots for MycRas and WT MEF data. 
  #
  # Parameters
  # ----------
  # XYZ_exp : string used to name first exp.
  # class_num : integer. number of classes to specify for functions with a k parameter.
  # method : string. method argument for hierarchical clustering. see hclust docs for options.
  # nb_clust : boolean. Run NbClust.
  # tsne_run : boolean. Run tsne and produce plots.
  #
  # Returns
  # -------
  # None.
  #
  
  # Load experiments from lists of experiment string locations
  myc_exp_sizes = c()
  myc_df <- data.frame()
  for (mexp in myc_exps){
    m_df = read.csv(paste('data/', mexp, '/', 'exp_motility_statistics.csv', sep =''))
    myc_exp_sizes = c(myc_exp_sizes, nrow(m_df))
    myc_df = rbind(myc_df, m_df)
  }
  
  wt_exp_sizes = c()
  wt_df <- data.frame()
  for (wtexp in wt_exps){
    wexp_df = read.csv(paste('data/', wtexp, '/', 'exp_motility_statistics.csv', sep =''))
    wt_exp_sizes = c(wt_exp_sizes, nrow(wexp_df))
    wt_df = rbind(wt_df, wexp_df)
  }
  all_exp_sizes <- c(myc_exp_sizes, wt_exp_sizes)
  df <- rbind(myc_df, wt_df)

  df <- df[ -c(1,2) ]
  df.raw <- df
  df <- data.frame( scale(df) )
  df.noturn <- df[,1:55]
  df.class <- c(rep("MycRas", nrow(myc_df)), rep('WT', nrow(wt_df)))
  df.class_num <- c(rep(1, nrow(myc_df)), rep(2, nrow(wt_df)))
  df.plate <- c()
  for (i in 1:length(all_exp_sizes)){
    df.plate <- c(df.plate, rep(i, all_exp_sizes[i]))
  }
  class_df <- df
  class_df$class <- as.factor(df.class)
  class_df.raw <- df.raw
  class_df.raw$class <- as.factor(df.class)
  
  if (sum(min_feat_set) > 1){
    df.noturn = df.noturn[,as.logical(as.matrix(min_feat_set))]
  } else if (sum(min_feat_set) == 1){
    f = colnames(df.noturn)[as.logical(as.matrix(min_feat_set))]
    df.noturn = data.frame(df.noturn[,as.logical(as.matrix(min_feat_set))])
    colnames(df.noturn) = f
  }
  
  stopifnot(nrow(df.noturn) == nrow(df.raw))
  
  df.pca_eigens <- prcomp(df.noturn)
  df.comp <- data.frame(df.pca_eigens$x[,1:min(c(30, ncol(df.noturn)))])
  
  if (ncol(df.noturn) == 1){
  } else {
    df.g <- heatmap3_plots(df.noturn, df.class, plot_path, class_num, method = method, prefix=prefix)
  }
  df.groups <- cutree(hclust(dist(df.comp[,1:min(c(30, ncol(df.noturn)))]), method = method), class_num)
  if (tsne_run){tsne_plots(df.noturn, df.groups, df.class, df.plate, plot_path, method = method, width = 2, height = 2, experiment = prefix)}

  df.cluster_membership <- mef_cluster_props(class_df, df.groups)
  mef_plot_cluster_memberships(df.cluster_membership, plot_path, prefix=prefix)
  mef_contingency_table_preftest(df.cluster_membership, df.groups, plot_path, experiment = paste(prefix, "MEF_ClustPref", sep=''))
  plot_cluster_features(df.noturn, df.groups, plot_path, anova_results = anova_results, width = 7, height = 4, prefix=prefix)

  capture.output(df.cluster_membership, file = paste(plot_path, prefix, "cluster_membership.txt", sep = ''))

  if (nb_clust == TRUE){
    clustval(df.comp, class_num = class_num, plot_path = plot_path, prefix = paste(prefix, "MEF_ClustVal_", sep=''), method=method)
  }
}

myc_exp = "mef/mycras/20160917"
myc2_exp = "mef/mycras/20160918"
wt_exp = "mef/wt/20160711"
wt2_exp = "mef/wt/20160925_0"
wt3_exp = "mef/wt/20160925_1"
wt4_exp = "mef/wt/20160927"
class_num = 3
method='ward.D2'
nb_clust = F
tsne_run = F

myc_exps <- c(myc_exp, myc2_exp)
wt_exps <- c(wt_exp, wt2_exp, wt3_exp, wt4_exp)

prcs = c(1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90)
for (i in 1:length(prcs)){
  p = prcs[i]
  min_feature_set = read.csv(paste('python/classif/mycwt_min_features_bool0', sprintf('%02d', p), '.csv', sep=''), header = F)
  mycras_wt_feature_reduction(myc_exps, wt_exps, class_num=class_num, method=method, 
                              nb_clust=nb_clust, tsne_run=tsne_run, 
                              prefix=paste('svm_min_feats', as.character(p), '_', sep=''), 
                              min_feat_set = min_feature_set)
}
