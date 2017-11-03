# MuSC : Myoblast comparisons

source('hm_common_fnx.R')
source('pseudotiming.R')

musc1 = 'musc/20160626'
musc2 = 'musc/20160701_0'
musc3 = 'musc/20160701_1'

myo1 <- "myoblast/20160623"
myo2 <- "myoblast/20160720"


load_muscs <- function(experiment1, experiment2, experiment3, method='ward.D', class_num=3){
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
  
  all_plates <- rbind(plate1, plate2, plate3)
  all_plates.class <- c(plate1.class, plate2.class, plate3.class)
  all_plates.class_num <- c(plate1.class_num, plate2.class_num, plate3.class_num)
  all_plates <- all_plates[ -c(1,2) ]
  all_plates.plate <- c( rep(1, nrow(plate1)), rep(2, nrow(plate2)), rep(3, nrow(plate3)))
  
  r <-  c(49, 84, 624, 1238, 1386, 1472, 1599, 1738, 1776, 1796, 1907, 1954, 1963,  2091, 2468, 2543, 2989, 3624, 3873, 4006, 4145, 4192)
  all_plates <- all_plates[-r,]
  all_plates.class <- all_plates.class[-r]
  all_plates.class_num <- all_plates.class_num[-r]
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
  
  all_plates.comp = prcomp(all_plates.noturn)$x
  all_plates.groups <- cutree(hclust(dist(all_plates.comp[,1:30], method='euclidean'), method=method), class_num)
  
  return(list(all_plates.raw, all_plates.class, all_plates.groups))
}

load_myoblasts <- function(myo_exp1, myo_exp2){
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

  return(list(df.raw, df.class))
}

plot_pca_group_vectors <- function(df.comp, df.groups, groupwise_vectors, plot_path, timestep='all', amp=7, width=4, height=4, exp=NULL){
  
  pca_theme <- theme(axis.text.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks=element_blank(),
                     legend.position="right",
                     panel.background=element_blank(),
                     panel.border=element_blank(),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),
                     plot.background=element_blank())
  
  nb_groups <- length(unique(df.groups))
  # Find centers of each group for start of arrows
  gcenters <- data.frame(matrix(nrow=nb_groups, ncol=2))
  colnames(gcenters) <- c('PC1', 'PC2')
  
  for (i in 1:nb_groups){
    g = unique(df.groups)[i]
    gpc = df.comp[df.groups == g,1:2]
    gcenters[i,] <- colMeans(gpc)
  }
  
  gvecs <- groupwise_vectors[groupwise_vectors$timestep == timestep,]
  gcenters$vPC1 <- gvecs$vdim00_mean
  gcenters$vPC2 <- gvecs$vdim01_mean*-1 # sklearn PCA inverts PC2 loadings
  gcenters$PC1end <- gcenters$PC1 + gcenters$vPC1
  gcenters$PC2end <- gcenters$PC2 + gcenters$vPC2
  
  
  p <- ggplot() + geom_point(data = df.comp, aes(x=PC1, y=PC2, color=df.groups)) 
  p <- p + geom_segment(data = gcenters, mapping=aes(x=PC1, xend=PC1+(vPC1*amp), y=PC2, yend=PC2+(vPC2*amp)), 
                        arrow = arrow(length = unit(0.15, 'cm'), type = 'closed'), size = 0.8)
  p <- p + pca_theme + scale_color_discrete('Cluster')
  ggsave(p, path = plot_path, filename = paste(exp, 'groupwise_transition_vectors_PCA.png', sep='_'), width = width, height = height)
  
}

plot_group_trans_mag <- function(groupwise_vectors, plot_path, timestep='all', ordering = c(1,3,2,4), exp=NULL, width=3, height = 3){

  
  gv_t <- groupwise_vectors[groupwise_vectors$timestep == timestep,]
  
  gv_t$mag = sqrt(gv_t$vdim00_mean^2 + gv_t$vdim01_mean^2)
  gv_t$mag_sem = sqrt(gv_t$vdim00_sem^2 + gv_t$vdim01_sem^2)
  gv_t$upper = gv_t$mag + gv_t$mag_sem
  gv_t$lower = gv_t$mag - gv_t$mag_sem
  gv_t$lower[gv_t$lower < 0] <- 0
  gv_t$ordering <- ordering
  
  #Plotting
  cluster_theme <- theme(
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank(),
    axis.line.x = element_line(size = 0.5),
    axis.line.y = element_line(size = 0.5),
    axis.text.x = element_text(size = 10), 
    text= element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )
  
  p <- ggplot() + geom_errorbar(data=gv_t, aes(x=ordering, ymin=lower, ymax=upper))
  p <- p + geom_bar(data=gv_t, aes(x=ordering, y=mag, fill=group), stat = 'identity') + scale_color_discrete('Cluster')
  p <- p + labs(x='Progression Order', y='Transition Magnitude (PC units)', title='State Specific Transition Rates')
  p <- p + cluster_theme + scale_y_continuous(expand = c(0,0), limits = c(0, 1.10*max(gv_t$upper)))
  
  ggsave(p, path=plot_path, filename = paste(exp, 'groupwise_trans_mag.png', sep='_'), width = width, height = height)
}


compare_musc_myo <- function(musc1, musc2, musc3, myo1, myo2, class_num, method='ward.D', nb_clust=F, tsne_run=F, pseudoT=F){
  plot_path = 'data/comparisons/muscmyo/'
  
  musc_list <- load_muscs(musc1, musc2, musc3)
  musc_df <- musc_list[[1]]
  musc_df.class <- musc_list[[2]]
  musc_df.groups <- musc_list[[3]]
  for (i in 1:max(musc_df.groups)){
    musc_df.groups[musc_df.groups==i] = paste('MuSC_', i, sep='')
  }
  
  myo_list <- load_myoblasts(myo1, myo2)
  myo_df <- myo_list[[1]]
  myo_df.class <- myo_list[[2]]
  
  df <- rbind(musc_df, myo_df)
  df.type <- c( rep('MuSC', nrow(musc_df)), rep('Myoblast', nrow(myo_df)) )
  df.raw <- df
  df <- data.frame(scale(df))
  df.noturn <- df[,1:55]
  
  df.groups <- c(musc_df.groups, rep('Myoblast', nrow(myo_df)))
  write.csv(data.frame(df.groups), file = paste(plot_path, 'muscmyo_hclust_labels.csv', sep=''), row.names = F)
  
  df.comp <- pca_plots(df.noturn, df.type, df.groups, plot_path, class_num)
  
  groupwise_vectors <- read.csv('data/comparisons/muscmyo/muscmyo_groupwise_transition_vectors.csv')
  plot_pca_group_vectors(df.comp, df.groups, groupwise_vectors, plot_path, timestep = 'all', amp=7, exp='MuSC_Myo')
  plot_group_trans_mag(groupwise_vectors, plot_path, exp = 'MuSC_Myo')
  
  if (tsne_run){tsne_plots(df.comp, df.groups, df.type, df.type, plot_path, perplexity = c(30,50,70), width=2.5, height=2.5)}
  
  
  if (pseudoT){
    ## Run Monocle
    # create CellDataSet using a subsample of musc cells to allow tractable computation
    idx <- sample(nrow(musc_df), nrow(myo_df)*4)
    #idx <- 1:nrow(musc_df) #uncomment to use all cells
    cds_df <- data.frame( scale( rbind( musc_df[idx,], myo_df ) ) )
    rownames(cds_df) <- 1:nrow(cds_df)
    cds_df.type <- c(df.type[idx], df.type[(nrow(musc_df)+1):length(df.type)])
    cds_df.groups <- c(df.groups[idx], df.groups[(nrow(musc_df)+1):length(df.groups)])
    cds <- prepare_cellDataSet(cds_df[,1:55], cds_df.type, cds_df.groups) # exclude turn statistics
    # Find ordering features as loadings of PC1
    df.pc <- princomp(df.noturn)
    aload <- abs( with(df.pc, unclass(loadings)) ) # get abs of all loadings
    pload <- sweep(aload, 2, colSums(aload), "/") # proportion per PC
    write.csv(pload, file = paste(plot_path, 'muscmyo_pca_loadings.csv', sep=''), row.names = F)
    ordering_features <- tail(sort(pload[,1]), 20) # top 20 loadings in PC1, get names with labels()
    cds <- mst_monocle(cds, ordering_features)
    plot_monocle(cds, plot_path, exp='MuSC_Myoblast_JointSpace')
  }

  
}


compare_musc_myo(musc1, musc2, musc3, myo1, myo2, class_num = 4, tsne_run = T, pseudoT = T)

