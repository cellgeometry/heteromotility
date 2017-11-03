# Local cell density calculation
require(ggplot2)
require(reshape2)
require(Rtsne)
source('hm_common_fnx.R')

## MuSC

class_num=3
method='ward.D'
plot_path='data/musc/lcd/'
experiment1 = 'musc/20160626'
experiment2 = 'musc/20160701_0'
experiment3 = 'musc/20160701_1'
class_num = 4
method = 'ward.D'
nb_clust = F
tsne_run = F

musc_lcd <- function(experiment1, experiment2, experiment3, class_num=4, method='ward.D'){
  plot_path = "data/musc/lcd/"
  
  fgf_df1 = read.csv( paste('data/', experiment1, '/', 'fgf2_exp_motility_statistics.csv', sep = '') )
  nofgf_df1 = read.csv( paste('data/', experiment1, '/', 'nofgf2_exp_motility_statistics.csv', sep = '') )
  fgf_dens1 = read.csv( paste('data/', experiment1, '/', 'fgf2_lcdensity.csv', sep = ''), header=F )
  nofgf_dens1 = read.csv( paste('data/', experiment1, '/', 'nofgf2_lcdensity.csv', sep = ''), header=F )
  
  fgf_df2 = read.csv( paste('data/', experiment2, '/', 'fgf2_exp_motility_statistics.csv', sep = '') )
  nofgf_df2 = read.csv( paste('data/', experiment2, '/', 'nofgf2_exp_motility_statistics.csv', sep = '') )
  fgf_dens2 = read.csv( paste('data/', experiment2, '/', 'fgf2_lcdensity.csv', sep = ''), header=F )
  nofgf_dens2 = read.csv( paste('data/', experiment2, '/', 'nofgf2_lcdensity.csv', sep = ''), header=F )
  
  fgf_df3 = read.csv( paste('data/', experiment3, '/', 'fgf2_exp_motility_statistics.csv', sep = '') )
  nofgf_df3 = read.csv( paste('data/', experiment3, '/', 'nofgf2_exp_motility_statistics.csv', sep = '') )
  fgf_dens3 = read.csv( paste('data/', experiment3, '/', 'fgf2_lcdensity.csv', sep = ''), header=F )
  nofgf_dens3 = read.csv( paste('data/', experiment3, '/', 'nofgf2_lcdensity.csv', sep = ''), header=F )
  
  plate1 = rbind(fgf_df1, nofgf_df1)
  plate1.class <- c(rep("FGF2+", nrow(fgf_df1)), rep("FGF2-", nrow(nofgf_df1)))
  plate1.class_num <- c(rep(1, nrow(fgf_df1)), rep(2, nrow(nofgf_df1)))
  plate1.dens <- rbind(fgf_dens1, nofgf_dens1)
  plate2 = rbind(fgf_df2, nofgf_df2)
  plate2.class <- c(rep("FGF2+", nrow(fgf_df2)), rep("FGF2-", nrow(nofgf_df2)))
  plate2.class_num <- c(rep(1, nrow(fgf_df2)), rep(2, nrow(nofgf_df2)))
  plate2.dens <- rbind(fgf_dens2, nofgf_dens2)
  plate3 = rbind(fgf_df3, nofgf_df3)
  plate3.class <- c(rep("FGF2+", nrow(fgf_df3)), rep("FGF2-", nrow(nofgf_df3)))
  plate3.class_num <- c(rep(1, nrow(fgf_df3)), rep(2, nrow(nofgf_df3)))
  plate3.dens <- rbind(fgf_dens3, nofgf_dens3)
  
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
  all_plates.dens <- rbind(plate1.dens, plate2.dens, plate3.dens)
  all_plates.dens <- all_plates.dens[-r]
  
  all_plates <- data.frame( scale(all_plates) )
  all_plates.noturn <- all_plates[,1:55]
  class_all_plates <- all_plates
  class_all_plates$class <- all_plates.class
  class_all_plates.noturn <- all_plates.noturn
  class_all_plates.noturn$class <- all_plates.class
  class_all_plates.raw <- all_plates.raw
  class_all_plates.raw$class <- all_plates.class
  
  # Cluster density into discrete bins
  density <- scale(all_plates.dens)
  #density[density == 0] <- NA
  density.k <- kmeans(na.omit(density), 3)
  #density.state <- rep(NA, length(density))
  density.state <- density.k$cluster
  
  dens_means <- data.frame(matrix(nrow=3, ncol=2))
  colnames(dens_means) <- c('Mean LCDI', 'Count')
  for (i in 1:3){
    dens_means[i,1] <- mean(density[density.state==i])
    dens_means[i,2] <- sum(density.state==i)
    
  }
  write.csv(dens_means, file = 'data/musc/lcd/density_clust_mean_vals.csv', row.names = F)
  
  all_plates.comp <- pca_plots(all_plates.noturn, all_plates.class, density.state, plot_path, class_num)
  all_plates.groups <- cutree(hclust(dist(all_plates.comp[,1:30], method='euclidean'), method=method), class_num)
  tsne_plots(all_plates.comp, all_plates.groups, all_plates.class, density.state, plot_path, perplexity = 50)
  
  density_props_df <- density_cluster_membership(all_plates.groups, density.state, plot_path = plot_path, prefix='MuSC_lcd')
  plot_density_cluster_membership(density_props_df, plot_path = plot_path, prefix = 'MuSC_lcd')
  
  density_reg_df <- density_regression_tests(all_plates.raw, density)
  write.csv(density_reg_df, file=paste(plot_path, 'MuSC_lcd', '_density_reg_tests.csv', sep=''), row.names = F)

}



## MEF
myc_exp = "mef/mycras/20160917"
myc2_exp = "mef/mycras/20160918"
wt_exp = "mef/wt/20160711"
wt2_exp = "mef/wt/20160925_0"
wt3_exp = "mef/wt/20160925_1"
wt4_exp = "mef/wt/20160927"
class_num = 3
method='ward.D2'

myc_exps <- c(myc_exp, myc2_exp)
wt_exps <- c(wt_exp, wt2_exp, wt3_exp, wt4_exp)

mycras_wt_density <- function(myc_exps, wt_exps, class_num, method='ward.D2', nb_clust=F, tsne_run=F, plot_path='data/mef/lcd/', prefix=''){
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
  plot_path = 'data/mef/lcd/'
  # Load experiments from lists of experiment string locations
  myc_exp_sizes = c()
  myc_df <- data.frame()
  myc_dens <- data.frame()
  for (mexp in myc_exps){
    m_df = read.csv(paste('data/', mexp, '/', 'exp_motility_statistics.csv', sep =''))
    m_dens <- read.csv(paste('data/', mexp, '/', 'exp_lcdensity.csv', sep =''), header = F)
    myc_exp_sizes = c(myc_exp_sizes, nrow(m_df))
    myc_df = rbind(myc_df, m_df)
    myc_dens <- rbind(myc_dens, m_dens)
  }
  
  wt_exp_sizes = c()
  wt_df <- data.frame()
  wt_dens <- data.frame()
  for (wtexp in wt_exps){
    wexp_df = read.csv(paste('data/', wtexp, '/', 'exp_motility_statistics.csv', sep =''))
    wexp_dens <- read.csv(paste('data/', wtexp, '/', 'exp_lcdensity.csv', sep =''), header = F)
    wt_exp_sizes = c(wt_exp_sizes, nrow(wexp_df))
    wt_df = rbind(wt_df, wexp_df)
    wt_dens = rbind(wt_dens, wexp_dens)
  }
  all_exp_sizes <- c(myc_exp_sizes, wt_exp_sizes)
  df <- rbind(myc_df, wt_df)
  density <- scale(rbind(myc_dens, wt_dens))
  
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
  
  density.k <- kmeans(na.omit(density), 3)
  density.state <- density.k$cluster
  
  dens_means <- data.frame(matrix(nrow=3, ncol=2))
  colnames(dens_means) <- c('Mean LCDI', 'Count')
  for (i in 1:3){
    dens_means[i,1] <- mean(density[density.state==i])
    dens_means[i,2] <- sum(density.state==i)
    
  }
  write.csv(dens_means, file = 'data/mef/lcd/density_clust_mean_vals.csv', row.names = F)
  
  df.comp <- pca_plots(df.noturn, df.class, density.state, plot_path, class_num, prefix=prefix)
  df.groups <- cutree(hclust(dist(df.comp[,1:30]), method = method), class_num)
  if (tsne_run){tsne_plots(df.noturn, df.groups, df.class, density.state, plot_path, method = method, width = 2, height = 2, experiment = prefix)}

  
  myc_idx <- df.class=='MycRas'
  wt_idx <- df.class=='WT'
  
  density_props_df <- density_cluster_membership(df.groups[myc_idx], density.state[myc_idx], plot_path = plot_path, prefix='MEF_MycRas_lcd')
  plot_density_cluster_membership(density_props_df, plot_path = plot_path, prefix = 'MEF_MycRas_lcd')
  
  density_reg_df <- density_regression_tests(df.raw[myc_idx,], density[myc_idx,])
  write.csv(density_reg_df, file=paste(plot_path, 'MEF_WT_lcd', '_density_reg_tests.csv', sep=''), row.names = F)  
  
  density_props_df <- density_cluster_membership(df.groups[wt_idx], density.state[wt_idx], plot_path = plot_path, prefix='MEF_WT_lcd')
  plot_density_cluster_membership(density_props_df, plot_path = plot_path, prefix = 'MEF_WT_lcd')
  
  density_reg_df <- density_regression_tests(df.raw[wt_idx,], density[wt_idx,])
  write.csv(density_reg_df, file=paste(plot_path, 'MEF_WT_lcd', '_density_reg_tests.csv', sep=''), row.names = F)  
}

## Myoblast

myo_exp1 <- "myoblast/20160623"
myo_exp2 <- "myoblast/20160720"
method = 'ward.D2'
class_num = 2
nb_clust = F
tsne_run = F

density_myoblasts <- function(myo_exp1, myo_exp2, plot_path, class_num, method='ward.D2', nb_clust=F, tsne_run=F){
  
  plot_path <- "data/myoblast/lcd/"
  
  fgf_df1 <- read.csv(paste('data/', myo_exp1, '/', 'fgf2_exp_motility_statistics.csv', sep = ''))
  nofgf_df1 <- read.csv(paste('data/', myo_exp1, '/', 'nofgf2_exp_motility_statistics.csv', sep = ''))
  fgf_dens1 <- read.csv(paste('data/', myo_exp1, '/', 'fgf2_lcdensity.csv', sep = ''), header = F)
  nofgf_dens1 <- read.csv(paste('data/', myo_exp1, '/', 'nofgf2_lcdensity.csv', sep = ''), header = F)
  
  fgf_df2 <- read.csv(paste('data/', myo_exp2, '/', 'fgf2_exp_motility_statistics.csv', sep = ''))
  nofgf_df2 <- read.csv(paste('data/', myo_exp2, '/', 'nofgf2_exp_motility_statistics.csv', sep = ''))
  fgf_dens2 <- read.csv(paste('data/', myo_exp2, '/', 'fgf2_lcdensity.csv', sep = ''), header = F)
  nofgf_dens2 <- read.csv(paste('data/', myo_exp2, '/', 'nofgf2_lcdensity.csv', sep = ''), header = F)
  
  plate1 = rbind(fgf_df1, nofgf_df1)
  plate1.class <- c(rep("FGF2+", nrow(fgf_df1)), rep("FGF2-", nrow(nofgf_df1)))
  plate1.class_num <- c(rep(1, nrow(fgf_df1)), rep(2, nrow(nofgf_df1)))
  plate1.dens <- rbind(fgf_dens1, nofgf_dens1)
  plate2 = rbind(fgf_df2, nofgf_df2)
  plate2.class <- c(rep("FGF2+", nrow(fgf_df2)), rep("FGF2-", nrow(nofgf_df2)))
  plate2.class_num <- c(rep(1, nrow(fgf_df2)), rep(2, nrow(nofgf_df2)))
  plate2.dens <- rbind(fgf_dens2, nofgf_dens2)
  
  r1 = c(6,15,19,43,48,49,51,53,66,70,115,135,150,154,167,202,204,205,208)
  r2 = c(3, 18, 19, 92, 97, 104, 112)
  
  plate1 = plate1[-r1,]
  plate1.class = plate1.class[-r1]
  plate1.class_num = plate1.class_num[-r1]
  plate1.dens <- plate1.dens[-r1,]
  plate2 = plate2[-r2,]
  plate2.class = plate2.class[-r2]
  plate2.class_num = plate2.class_num[-r2]
  plate2.dens <- plate2.dens[-r2,]
  
  
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
  
  density <- c(plate1.dens, plate2.dens)
  density.k <- kmeans(na.omit(density), 3)
  density.state <- density.k$cluster
  
  df.comp <- pca_plots(df.noturn, df.class, density.state, plot_path, class_num)
  df.groups <- cutree(hclust(dist(df.comp[,1:30]), method=method), class_num)
  if (tsne_run){tsne_plots(df.comp, df.groups, df.class, density.state, plot_path, width=2.5, height=2.5)}

  density_props_df <- density_cluster_membership(df.groups, density.state, plot_path = plot_path, prefix='Myoblast_lcd')
  plot_density_cluster_membership(density_props_df, plot_path = plot_path, prefix = 'Myoblast_lcd')
  
  density_reg_df <- density_regression_tests(df.raw, density[,1])
  write.csv(density_reg_df, file=paste(plot_path, 'Myoblast_lcd', '_density_reg_tests.csv', sep=''), row.names = F)  
  rsq_hist <- ggplot(density_reg_df, aes(Rsq)) + geom_histogram(bins=100) + labs(title='Feature ~ Density Relationships', x='R^2', y='Frequency')
}


process_myoblasts(myo_exp1, myo_exp2, plot_path, class_num=class_num, method=method, nb_clust=nb_clust)


plot_all_rsq <- function(){
  myoblast <- read.csv('data/myoblast/lcd/Myoblast_lcd_density_reg_tests.csv')
  mycras <- read.csv('data/mef/lcd/MEF_MycRas_lcd_density_reg_tests.csv')
  wt <- read.csv('data/mef/lcd/MEF_WT_lcd_density_reg_tests.csv')
  musc <- read.csv('data/musc/lcd/MuSC_lcd_density_reg_tests.csv')
  
  df <- rbind(myoblast[1:55,], mycras[1:55,], wt[1:55,], musc[1:55,])
  df$CellType <- c( rep('Myoblast', 55), rep('MEF MycRas', 55), rep('MEF WT', 55), rep('MuSC', 55))
  
  th <- theme(
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
  
  p <- ggplot(df, aes(Rsq, fill=CellType)) + geom_histogram(bins=100)
  p <- p + labs(title='Feature ~ Density Relationships', x='R-squared Value', y='Number of Features')
  p <- p + scale_y_continuous(expand = c(0,0), limits = c(0, 50)) + th
  ggsave(p, filename='data/comparisons/local_cell_density_rsq.png', width=4, height=3, units='in')
}
