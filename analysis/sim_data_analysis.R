levy_df <- read.csv('data/sims/sim_stats/motility_statistics_pf_t100_mu5.csv')
brw_df <- read.csv('data/sims/sim_stats/motility_statistics_brw_t100_mu5.csv')
rw_df <- read.csv('data/sims/sim_stats/motility_statistics_rw_t100_mu5.csv')
fbm_df <- read.csv('data/sims/sim_stats/motility_statistics_fbm_t100_mu5.csv')
plot_path <- "data/sims/sim_plots/"

levy_dfs <- levy_df[2000:3000,]
brw_dfs <- brw_df[2000:3000,]
rw_dfs <- rw_df[2000:3000,]
fbm_dfs <- fbm_df[2000:3000,]

df <- rbind(levy_dfs, brw_dfs, rw_dfs, fbm_dfs)
df <- df[ -c(1,2) ]
df <- data.frame(scale(df))
df.class <- c( rep("Levy", nrow(levy_dfs)), rep("Biased RW", nrow(brw_dfs)), rep("Unbiased RW", nrow(rw_dfs)), rep("fBm", nrow(fbm_dfs)) )
df.class_num <- c( rep(1, nrow(levy_dfs)), rep(2, nrow(brw_dfs)), rep(3, nrow(rw_dfs)), rep(4, nrow(fbm_dfs)) )

# PCA loadings
df.pc <- princomp(df)
df.pc.loadings <- loadings(df.pc)
df.pc.summary <- summary(df.pc)

# PCA plot of first 3 PCs
# use k-means clusters for colors
df.pca_eigens <- prcomp(df)
df.comp <- data.frame(df.pca_eigens$x[,1:30])
df.k <- kmeans(df.comp, 4, nstart=25, iter.max=1000)

# PCA Plot
require(ggplot2)
no_bg <- theme(#axis.line.x = element_line(color="black", size = 0.1),
  #axis.line.y = element_line(color="black", size = 0.1),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  legend.position="right",
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  plot.background=element_blank())

pca_class <- ggplot(df.comp, aes(x = PC1, y = PC2,  col = as.factor(df.class))) + geom_point() + no_bg + scale_colour_discrete(name = "Model Class")
ggsave(pca_class, filename = "pca_class.png", path = plot_path, width = 5, height = 5, units = "in")

# Heirarchical clustering
df.noturn <- df[,1:56]
df.noturn.mt <- as.matrix(df.noturn)

df.d <- dist(df.noturn, method="euclidean")
df.fit <- hclust(df.d, method="ward.D")
df.groups <- cutree(df.fit, k = 4) # 4 actual ground truth classes

capture.output(table(df.groups,df.class), file = paste(plot_path, "hclust_accuracy.txt", sep =''))

library(gplots)


hclust_wardD <- function(x) hclust(x, method='ward.D')

png(filename = paste(plot_path, "heatmap_wardD_classlab.png", sep=''), height = 10, width = 8, units = "in", res = 300)
heatmap.2(df.noturn.mt, main = "Simulated Motility Models", hclustfun = hclust_wardD, col=redgreen, scale = 'row', trace = "none", dendrogram = "row", density.info = "none", RowSideColors = as.character(df.class_num), labRow = NULL, margins = c(10,10))
dev.off()

# TSNE

require(Rtsne)
df.noturn.tsne <- Rtsne(as.matrix(df.comp), perplexity = 30)

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

t_class <- ggplot(data.frame(df.noturn.tsne$Y), aes(x=X1, y=X2, col = as.factor(df.class))) + geom_point(size=1) + labs(title = "Simulated Motility Models t-SNE") + no_axes_bg
ggsave(t_class, path = plot_path, filename = "tsne_simulations.png", width = 5, height = 5, units = "in")

# ICA
require(ica)
df.noturn.ica <- icafast(df.noturn, nc =4, maxit= 10000, alpha = 1.5)
plot(df.noturn.ica$S, col = as.factor(df.class))
df.ica_points <- data.frame(df.noturn.ica$S)
colnames(df.ica_points) <- c("IC1", "IC2", "IC3", "IC4")
ica_class <- ggplot(df.ica_points, aes(x = IC1, y = IC2, col = as.factor(df.class))) + geom_point(size=1) + no_bg + scale_colour_discrete(name = "Model Class")
ggsave(ica_class, path = plot_path, filename = "ica_class.png", width = 5, height = 5, units = "in")

# Cluster Validation

library(cluster)
df.si <- silhouette(df.class_num, dist = df.d)
plot(df.si)

mv <- manova(as.matrix(df.comp) ~ as.factor(df.class))
capture.output(summary(mv), file = paste(plot_path, "manova_pca_class.txt", sep=''))


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
capture.output(meanAcc, file = paste(plot_path, 'RFclassification_5CV_acc.txt', sep=''))
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

conf_df$Count = as.numeric(conf_df$Count)

heat_theme <- theme(
                   axis.ticks=element_blank(),
                   legend.position="right",
                   panel.background=element_blank(),
                   panel.border=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),
                   plot.background=element_blank())

p <- ggplot(data=conf_df, aes(x=Predicted_Class, y=True_Class)) + geom_tile(aes(fill=Count)) 
p <- p + scale_fill_gradient() + geom_text(aes(label = sprintf("%1.0f", Count), colour='w'), vjust = 1) + scale_color_manual(values = c(w='white'))
p <- p + heat_theme + coord_fixed()

ggsave(p, filename='classification_conf_mat.png', path=plot_path, width=5, height=4)

# Save feature importance 

m <- randomForest(x=df.noturn.mt, y=as.factor(df.class), ntree=100, importance=T)
m.importance <- data.frame(m$importance)
sorted_importance <- m.importance[order(m.importance$MeanDecreaseAccuracy, decreasing = T),]
tab <- data.frame(cbind(rownames(sorted_importance)[1:10], sorted_importance[1:10,5]))
colnames(tab) <- c('Feature', 'Decrease in Accuracy')
write.csv(tab, file = paste(plot_path, 'rf_feature_importance.csv', sep=''), row.names = F)

