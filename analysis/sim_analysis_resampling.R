# Simulated data analysis with resampling and varied track length
source('hm_common_fnx.R')

plot_path <- "data/sims/sim_plots/"
class_num = 4
ns = c(100,500,1000,5000) # sample sizes
ts = c(50,100,500) # time lengths

no_bg <- theme(axis.line=element_blank(),
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

sim_analysis <- function(n=1000, t=100, plot_path='data/sims/sim_plots', class_num=4, prefix=''){
  # Loads a data.frame of motility statistics from simulated paths and performs analysis based
  # on a random sample of size `n` with track_length `t`
  # Performs tSNE, PCA, ICA, RF Classification
  
  # Load df's
  levy_df <- read.csv(paste('data/sims/sim_stats/motility_statistics_', 'pf_t', sprintf('%03d', t), '_mu5.csv', sep=''))
  brw_df <- read.csv(paste('data/sims/sim_stats/motility_statistics_', 'brw_t', sprintf('%03d', t), '_mu5.csv', sep=''))
  rw_df <- read.csv(paste('data/sims/sim_stats/motility_statistics_', 'rw_t', sprintf('%03d', t), '_mu5.csv', sep=''))
  fbm_df <- read.csv(paste('data/sims/sim_stats/motility_statistics_', 'fbm_t', sprintf('%03d', t), '_mu5.csv', sep=''))
  
  # sample `n` objects from each df
  levy_dfs <- levy_df[sample(1:nrow(levy_df), size=n, replace=FALSE),]
  brw_dfs <- brw_df[sample(1:nrow(brw_df), size=n, replace=FALSE),]
  rw_dfs <- rw_df[sample(1:nrow(rw_df), size=n, replace=FALSE),]
  fbm_dfs <- fbm_df[sample(1:nrow(fbm_df), size=n, replace=FALSE),]
  
  # build common df
  df <- rbind(levy_dfs, brw_dfs, rw_dfs, fbm_dfs)
  df <- df[ -c(1,2) ]
  df <- data.frame(scale(df))
  df.class <- c( rep("Levy", nrow(levy_dfs)), rep("Biased RW", nrow(brw_dfs)), rep("Unbiased RW", nrow(rw_dfs)), rep("fBm", nrow(fbm_dfs)) )
  df.class_num <- c( rep(1, nrow(levy_dfs)), rep(2, nrow(brw_dfs)), rep(3, nrow(rw_dfs)), rep(4, nrow(fbm_dfs)) )
  df.noturn <- df[,1:56]
  
  # dim reduction plots
  df.comp <- pca_plots(df.noturn, df.class, df.class, plot_path, class_num = class_num, clust_pcs=10, prefix=paste(prefix, 'n', sprintf('%05d', n), '_t', sprintf('%03d', t), '_', sep=''))

  # PCA
  df.pca_eigens <- prcomp(df.noturn)
  df.comp <- data.frame(df.pca_eigens$x[,1:30])
  df.pc <- princomp(df.noturn)
  aload <- abs( with(df.pc, unclass(loadings)) ) # get abs of all loadings
  pload <- sweep(aload, 2, colSums(aload), "/") # proportion per PC
  write.csv(pload, file = paste(plot_path, prefix, 'n', sprintf('%05d', n), '_t', sprintf('%03d', t), '_pca_loadings.csv', sep=''), row.names = F)

  # hierarchical clustering
  df.noturn.mt <- as.matrix(df.noturn)
  df.d <- dist(df.comp, method="euclidean")
  df.fit <- hclust(df.d, method="ward.D")

  df.groups <- cutree(df.fit, k = class_num) # 4 actual ground truth classes
  hclust_table <- table(df.groups,df.class)
  hclust_accuracy <- sum(apply(FUN=max, hclust_table, 2)/n)/ncol(hclust_table)
  capture.output(hclust_table, file = paste(plot_path, prefix, 'n', sprintf('%05d', n), '_t', sprintf('%03d', t), "_hclust_confusion.txt", sep =''))
  capture.output(hclust_accuracy, file = paste(plot_path, prefix, 'n', sprintf('%05d', n), '_t', sprintf('%03d', t), "_hclust_accuracy.txt", sep =''))
  
  # tSNE
  tsne_plots(df.comp, df.groups, df.class, df.class, plot_path, experiment = paste(prefix, 'n', sprintf('%05d', n), '_t', sprintf('%03d', t), '_', sep=''))
  
  
  library(gplots)
  hclust_wardD <- function(x) hclust(x, method='ward.D')
  png(filename = paste(plot_path, prefix, 'n', sprintf('%05d', n), '_t', sprintf('%03d', t), "_heatmap_wardD_classlab.png", sep=''), height = 10, width = 8, units = "in", res = 300)
  heatmap.2(df.noturn.mt, main = "Simulated Motility Models", hclustfun = hclust_wardD, col=redgreen, scale = 'row', trace = "none", dendrogram = "row", density.info = "none", RowSideColors = as.character(df.class_num), labRow = NULL, margins = c(10,10))
  dev.off()
  
  require(ica)
  df.noturn.ica <- icafast(df.noturn, nc =4, maxit= 1000, alpha = 1.5)
  plot(df.noturn.ica$S, col = as.factor(df.class))
  df.ica_points <- data.frame(df.noturn.ica$S)
  colnames(df.ica_points) <- c("IC1", "IC2", "IC3", "IC4")
  ica_class <- ggplot(df.ica_points, aes(x = IC1, y = IC2, col = as.factor(df.class))) + geom_point(size=1) + no_bg + scale_colour_discrete(name = "Model Class")
  ggsave(ica_class, path = plot_path, filename = paste(prefix, 'n', sprintf('%05d', n), '_t', sprintf('%03d', t), "_ica_class.png", sep=''), width = 5, height = 5, units = "in")
  
  # Cluster Validation
  library(cluster)
  df.si <- silhouette(df.class_num, dist = df.d)
  plot(df.si)
  

  mv <- manova(as.matrix(df.comp) ~ as.factor(df.class))
  capture.output(summary(mv), file = paste(plot_path, prefix, 'n', sprintf('%05d', n), '_t', sprintf('%03d', t), "_manova_pca_class.txt", sep=''))
  
  
  # RF classification
  require(randomForest)
  require(caret)
  
  cvRF <- function(x, y, k, ntree=100){
    folds <- createFolds(as.factor(y), k=k, list=T)
    nb_class <- length(unique(y))
    m_list <- c()
    acc <- c()
    conf_mat <- matrix(ncol=nb_class+1, nrow=nb_class, 0)
    
    for (i in 1:k){
      print(paste('Training Fold ', i, sep=''))
      x_test <- x[folds[[i]], ]
      y_test <- y[folds[[i]]]
      
      x_train <- x[-(folds[[i]]), ]
      y_train <- y[-(folds[[i]])]
      
      
      m <- randomForest(x=x_train, y=as.factor(y_train), ntree = ntree)
      m_list <- c(m_list, m)
      
      yhat <- predict(m, x_test)
      acc <- c(acc, (sum(yhat==y_test)/length(y_test)))
      conf_mat <- conf_mat + as.matrix(m$confusion)
      
    }
    
    conf_mat <- conf_mat/k
    return(list(conf_mat, acc))
  }
  
  cvList <- cvRF(df.noturn.mt, df.class, k=5, ntree=100)
  meanAcc <- mean(cvList[[2]])
  capture.output(meanAcc, file = paste(plot_path, prefix, 'n', sprintf('%05d', n), '_t', sprintf('%03d', t),'_RFclassification_5CV_acc.txt', sep=''))
  conf_mat <- cvList[[1]]
  conf_df <- data.frame(matrix(nrow=nrow(conf_mat)**2, ncol=3))
  colnames(conf_df) <- c('True_Class', 'Predicted_Class', 'Count')
  conf_df$Count <- as.numeric(conf_df$Count)
  
  k = 1
  for (i in 1:nrow(conf_mat)){
    
    for (j in 1:nrow(conf_mat)){
      conf_df[k,] <- c(rownames(conf_mat)[i], colnames(conf_mat)[j], conf_mat[i,j])
      k <- k + 1  
    }
    
  }
  
  heat_theme <- theme(
    axis.ticks=element_blank(),
    legend.position="right",
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank())
  
  conf_df$Count <- as.numeric(conf_df$Count)
  
  p <- ggplot(data=conf_df, aes(x=Predicted_Class, y=True_Class)) + geom_tile(aes(fill=Count)) 
  p <- p + scale_fill_gradient() + geom_text(aes(label = sprintf("%1.0f", Count), colour='w'), vjust = 1) + scale_color_manual(values = c(w='white'))
  p <- p + heat_theme + coord_fixed()
  
  ggsave(p, filename=paste(prefix, 'n', sprintf('%05d', n), '_t', sprintf('%03d', t), '_classification_conf_mat.png', sep=''), path=plot_path, width=5, height=4)
  
  # Save feature importance 
  m <- randomForest(x=df.noturn.mt, y=as.factor(df.class), ntree=100, importance=T)
  m.importance <- data.frame(m$importance)
  sorted_importance <- m.importance[order(m.importance$MeanDecreaseAccuracy, decreasing = T),]
  tab <- data.frame(cbind(rownames(sorted_importance)[1:10], sorted_importance[1:10,5]))
  colnames(tab) <- c('Feature', 'Decrease in Accuracy')
  write.csv(tab, file = paste(plot_path, prefix, 'n', sprintf('%05d', n), '_t', sprintf('%03d', t), '_rf_feature_importance.csv', sep=''), row.names = F)
  
}

# Analyze across sample sizes and time scales
for (i in 1:length(ns)){
  for (j in 1:length(ts)){
    for (k in 1:3){
      sim_analysis(n=ns[i], t=ts[j], plot_path = plot_path, class_num=4, prefix=paste('resamp_', sprintf('%02d', k), '_', sep=''))
    }
  }
}

# take average accuracy for each resampling run
acc_df <- read.csv('data/sims/sim_plots/resamp_aggregate_hclust_accuracy.csv')

plot_acc <- function(acc_df, type='hclust'){
  gg_blue <- '#00BFC4'
  
  sem <- function(x){sd(x)/sqrt(length(x))}
  t_acc_df <- acc_df[-c(2)]
  n_acc_df <- acc_df[-c(3)]
  
  t_means <- aggregate(. ~ t, data = t_acc_df, FUN = mean)
  t_stderr <- aggregate(. ~ t, data = t_acc_df, FUN = sem)
  t_means.m <- melt(t_means, value.name = "Mean", id.vars=c('t'))
  t_stderr.m <- melt(t_stderr, value.name = "SE", id.vars=c('t'))
  t_plt_df <- merge(t_means.m, t_stderr.m)
  t_plt_df$t <- as.factor(as.matrix(t_plt_df$t))
  
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
  
  t_mean_limits <- aes(ymax = t_plt_df$Mean + t_plt_df$SE, ymin = t_plt_df$Mean - t_plt_df$SE)
  t_mean <- ggplot(t_plt_df, aes(t, Mean, fill=t)) + geom_bar(position = "dodge", stat = "identity") 
  t_mean <- t_mean + geom_errorbar(t_mean_limits, position = position_dodge(0.9), width = 0.5)
  t_mean <- t_mean + scale_y_continuous(expand=c(0,0), limits = c(0,1)) + cluster_theme
  t_mean <- t_mean + labs(title = "Classification Accuracy", x="Path length [t]", y = "Accuracy (+/- SE)")
  ggsave(t_mean, filename = paste('data/sims/sim_plots/resamp_', type, '_accuracy_t.png', sep=''), width = 4, height = 3, units = 'in')
  
  n_means <- aggregate(. ~ n, data = n_acc_df, FUN = mean)
  n_stderr <- aggregate(. ~ n, data = n_acc_df, FUN = sem)
  n_means.m <- melt(n_means, value.name = "Mean", id.vars=c('n'))
  n_stderr.m <- melt(n_stderr, value.name = "SE", id.vars=c('n'))
  n_plt_df <- merge(n_means.m, n_stderr.m)
  n_plt_df$n <- as.factor(as.matrix(n_plt_df$n))
  
  n_mean_limits <- aes(ymax = n_plt_df$Mean + n_plt_df$SE, ymin = n_plt_df$Mean - n_plt_df$SE)
  n_mean <- ggplot(n_plt_df, aes(n, Mean, fill=n)) + geom_bar(position = "dodge", stat = "identity") 
  n_mean <- n_mean + geom_errorbar(n_mean_limits, position = position_dodge(0.9), width = 0.5)
  n_mean <- n_mean + scale_y_continuous(expand = c(0,0), limits = c(0,1)) + cluster_theme
  n_mean <- n_mean + labs(title = "Classification Accuracy", x="Sample size", y = "Accuracy (+/- SE)")
  ggsave(n_mean, filename = paste('data/sims/sim_plots/resamp_', type, '_accuracy_n.png', sep=''), width = 4, height = 3, units = 'in')
}

plot_acc(acc_df, type='hclust')
cv_acc_df <- read.csv('data/sims/sim_plots/resamp_aggregate_RF_CV_accuracy.csv')
plot_acc(cv_acc_df, type='RFCV')

# Analyze different parameter inputs
brws <- dir('data/sims/sim_stats/var_stats/', pattern='*brw*')
fbms <- dir('data/sims/sim_stats/var_stats/', pattern='*fbm*')
pfs <- dir('data/sims/sim_stats/var_stats/', pattern='*pf*')
urws <- dir('data/sims/sim_stats/var_stats/', pattern='*_rw*')

sim_var_analysis <- function(brw, fbm, pf, urw, n=1000, plot_path='data/sims/sim_plots/var_plots/', class_num=4, prefix=''){
  # Loads a data.frame of motility statistics from simulated paths and performs analysis based
  # on a random sample of size `n` with track_length `t`
  # Performs tSNE, PCA, ICA, RF Classification
  
  # Load df's
  brw_df <- read.csv(paste('data/sims/sim_stats/var_stats/', brw, sep=''))
  fbm_df <- read.csv(paste('data/sims/sim_stats/var_stats/', fbm, sep=''))
  levy_df <- read.csv(paste('data/sims/sim_stats/var_stats/', pf, sep=''))
  rw_df <- read.csv(paste('data/sims/sim_stats/var_stats/', urw, sep=''))
  
  # sample `n` objects from each df
  levy_dfs <- levy_df[sample(1:nrow(levy_df), size=n, replace=FALSE),]
  brw_dfs <- brw_df[sample(1:nrow(brw_df), size=n, replace=FALSE),]
  rw_dfs <- rw_df[sample(1:nrow(rw_df), size=n, replace=FALSE),]
  fbm_dfs <- fbm_df[sample(1:nrow(fbm_df), size=n, replace=FALSE),]
  
  # build common df
  df <- rbind(levy_dfs, brw_dfs, rw_dfs, fbm_dfs)
  df <- df[ -c(1,2) ]
  df <- data.frame(scale(df))
  df.class <- c( rep("Levy", nrow(levy_dfs)), rep("Biased RW", nrow(brw_dfs)), rep("Unbiased RW", nrow(rw_dfs)), rep("fBm", nrow(fbm_dfs)) )
  df.class_num <- c( rep(1, nrow(levy_dfs)), rep(2, nrow(brw_dfs)), rep(3, nrow(rw_dfs)), rep(4, nrow(fbm_dfs)) )
  df.noturn <- df[,1:56]
  
  # dim reduction plots
  df.comp <- pca_plots(df.noturn, df.class, df.class, plot_path, class_num = class_num, clust_pcs=10, prefix=paste(prefix, sep=''))
  
  # PCA
  df.pca_eigens <- prcomp(df.noturn)
  df.comp <- data.frame(df.pca_eigens$x[,1:30])
  df.pc <- princomp(df.noturn)
  aload <- abs( with(df.pc, unclass(loadings)) ) # get abs of all loadings
  pload <- sweep(aload, 2, colSums(aload), "/") # proportion per PC
  write.csv(pload, file = paste(plot_path, prefix, '_pca_loadings.csv', sep=''), row.names = F)
  
  # hierarchical clustering
  df.noturn.mt <- as.matrix(df.noturn)
  df.d <- dist(df.comp, method="euclidean")
  df.fit <- hclust(df.d, method="ward.D")
  
  df.groups <- cutree(df.fit, k = class_num) # 4 actual ground truth classes
  hclust_table <- table(df.groups,df.class)
  hclust_accuracy <- sum(apply(FUN=max, hclust_table, 2)/n)/ncol(hclust_table)
  capture.output(hclust_table, file = paste(plot_path, prefix, "_hclust_confusion.txt", sep =''))
  capture.output(hclust_accuracy, file = paste(plot_path, prefix, "_hclust_accuracy.txt", sep =''))
  
  # tSNE
  tsne_plots(df.comp, df.groups, df.class, df.class, plot_path, experiment = paste(prefix, sep=''))
  
  
  library(gplots)
  hclust_wardD <- function(x) hclust(x, method='ward.D')
  png(filename = paste(plot_path, prefix, "_heatmap_wardD_classlab.png", sep=''), height = 10, width = 8, units = "in", res = 300)
  heatmap.2(df.noturn.mt, main = "Simulated Motility Models", hclustfun = hclust_wardD, col=redgreen, scale = 'row', trace = "none", dendrogram = "row", density.info = "none", RowSideColors = as.character(df.class_num), labRow = NULL, margins = c(10,10))
  dev.off()
  
  require(ica)
  df.noturn.ica <- icafast(df.noturn, nc =4, maxit= 1000, alpha = 1.5)
  plot(df.noturn.ica$S, col = as.factor(df.class))
  df.ica_points <- data.frame(df.noturn.ica$S)
  colnames(df.ica_points) <- c("IC1", "IC2", "IC3", "IC4")
  ica_class <- ggplot(df.ica_points, aes(x = IC1, y = IC2, col = as.factor(df.class))) + geom_point(size=1) + no_bg + scale_colour_discrete(name = "Model Class")
  ggsave(ica_class, path = plot_path, filename = paste(prefix, "ica_class.png", sep=''), width = 5, height = 5, units = "in")
  
  # Cluster Validation
  library(cluster)
  df.si <- silhouette(df.class_num, dist = df.d)
  plot(df.si)
  
  
  mv <- manova(as.matrix(df.comp) ~ as.factor(df.class))
  capture.output(summary(mv), file = paste(plot_path, prefix, "manova_pca_class.txt", sep=''))
  
  
  # RF classification
  require(randomForest)
  require(caret)
  
  cvRF <- function(x, y, k, ntree=100){
    folds <- createFolds(as.factor(y), k=k, list=T)
    nb_class <- length(unique(y))
    m_list <- c()
    acc <- c()
    conf_mat <- matrix(ncol=nb_class+1, nrow=nb_class, 0)
    
    for (i in 1:k){
      print(paste('Training Fold ', i, sep=''))
      x_test <- x[folds[[i]], ]
      y_test <- y[folds[[i]]]
      
      x_train <- x[-(folds[[i]]), ]
      y_train <- y[-(folds[[i]])]
      
      
      m <- randomForest(x=x_train, y=as.factor(y_train), ntree = ntree)
      m_list <- c(m_list, m)
      
      yhat <- predict(m, x_test)
      acc <- c(acc, (sum(yhat==y_test)/length(y_test)))
      conf_mat <- conf_mat + as.matrix(m$confusion)
      
    }
    
    conf_mat <- conf_mat/k
    return(list(conf_mat, acc))
  }
  
  cvList <- cvRF(df.noturn.mt, df.class, k=5, ntree=100)
  meanAcc <- mean(cvList[[2]])
  capture.output(meanAcc, file = paste(plot_path, prefix, 'RFclassification_5CV_acc.txt', sep=''))
  conf_mat <- cvList[[1]]
  conf_df <- data.frame(matrix(nrow=nrow(conf_mat)**2, ncol=3))
  colnames(conf_df) <- c('True_Class', 'Predicted_Class', 'Count')
  conf_df$Count <- as.numeric(conf_df$Count)
  
  k = 1
  for (i in 1:nrow(conf_mat)){
    
    for (j in 1:nrow(conf_mat)){
      conf_df[k,] <- c(rownames(conf_mat)[i], colnames(conf_mat)[j], conf_mat[i,j])
      k <- k + 1  
    }
    
  }
  
  heat_theme <- theme(
    axis.ticks=element_blank(),
    legend.position="right",
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank())
  
  conf_df$Count <- as.numeric(conf_df$Count)
  
  p <- ggplot(data=conf_df, aes(x=Predicted_Class, y=True_Class)) + geom_tile(aes(fill=Count)) 
  p <- p + scale_fill_gradient() + geom_text(aes(label = sprintf("%1.0f", Count), colour='w'), vjust = 1) + scale_color_manual(values = c(w='white'))
  p <- p + heat_theme + coord_fixed()
  
  ggsave(p, filename=paste(prefix, 'classification_conf_mat.png', sep=''), path=plot_path, width=5, height=4)
  
  # Save feature importance 
  m <- randomForest(x=df.noturn.mt, y=as.factor(df.class), ntree=100, importance=T)
  m.importance <- data.frame(m$importance)
  sorted_importance <- m.importance[order(m.importance$MeanDecreaseAccuracy, decreasing = T),]
  tab <- data.frame(cbind(rownames(sorted_importance)[1:10], sorted_importance[1:10,5]))
  colnames(tab) <- c('Feature', 'Decrease in Accuracy')
  write.csv(tab, file = paste(plot_path, prefix, 'rf_feature_importance.csv', sep=''), row.names = F)
  
}

# Randomly select combinations of simulations with parameters varied

epochs = 10
var_groups_log <- data.frame(matrix(nrow=epochs, ncol=4))
colnames(var_groups_log) <- c('brw', 'fbm', 'pf', 'urw')

for (i in 1:epochs){
  
  brw <- sample(brws, 1)
  fbm <- sample(fbms, 1)
  pf <- sample(pfs, 1)
  urw <- sample(urws, 1)
  
  var_groups_log[i,] <- c(brw, fbm, pf, urw)
  sim_var_analysis(brw, fbm, pf, urw, prefix=paste('var', sprintf('%02g', i), '_', sep=''))
  
}

write.csv(var_groups_log, paste('data/sims/sim_stats/var_stats/', '00_var_groups_log.csv', sep=''), row.names=F)

# plot var mean acc 

var_acc_df_uc <- read.csv('data/sims/sim_plots/var_plots/var_aggregate_hclust_accuracy.csv')
var_acc_df_rf <- read.csv('data/sims/sim_plots/var_plots/var_aggregate_RF_CV_accuracy.csv')

cluster_theme <- theme(
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  plot.background=element_blank(),
  axis.line.x = element_line(size = 0.5),
  axis.line.y = element_line(size = 0.5),
  axis.text.x = element_text(size = 10, angle=30, hjust = 1), 
  text= element_text(size = 10)
)

plt_df <- data.frame(matrix(ncol=3, nrow=2))
plt_df[1,] <- c(mean(var_acc_df_uc$acc), sem(var_acc_df_uc$acc), 'Unsupervised Clustering')
plt_df[2,] <- c(mean(var_acc_df_rf$acc), sem(var_acc_df_rf$acc), 'Random Forest')
colnames(plt_df) <- c('Mean', 'SEM', 'Classification')
plt_df$Mean <- as.numeric(plt_df$Mean)
plt_df$SEM <- as.numeric(plt_df$SEM)


p <- ggplot(plt_df, aes(x=Classification, y=Mean, fill=Classification)) + geom_bar(stat='identity')
p <- p + geom_errorbar(aes(ymin=Mean-SEM, ymax=Mean+SEM, width=0.4))
p <- p + cluster_theme + scale_y_continuous(expand=c(0,0), limits = c(0, 1.1))
p <- p + labs(x='Classification Type', y='Accuracy (+/- SE)')

ggsave(p, path=plot_path, filename='var_unsup_rf_accuracies.png', width=4.5, height=3, units='in')
