## Predefined Cluster PFA 

predef_state_locations <- function(df, class_num, method='ward.D'){
  # Generates a matrix of state locations based on hierarchical clustering
  # of subpaths in df
  #
  # Parameters
  # ----------
  # df : data.frame
  #   N x M feature matrix, where cols 1:3 are IDs and 4:55 are features
  #   column 'cell_id' is formated 'C-T' where C and T are integers, C 
  #   specifying a unique cell_id while T is the timestep for the subpath
  # class_num : integer.
  #   number of classes for clustering.
  # method : string.
  #   method for hierarchical clustering using hclust()
  #
  # Returns
  # -------
  # predef_sl : data.frame
  #   N*T x 4 array, plate, Well.Xy, cell_id$cell, cell_id$split
  #   where N is the number of cells, T is the number of timesteps
  require(reshape2)
  df.vals <- data.frame( scale(df[,4:55]) )
  df.d <- dist(df.vals)
  df.fit <- hclust(df.d, method = method)
  df.groups <- cutree(df.fit, k = class_num)
  
  group_frame <- data.frame(df[,1:3], df.groups)
  group_frame <- transform(group_frame, cell_id = colsplit(cell_id, "-", names = c("cell", "split")))
  
  return(group_frame)
}

sl_contingency <- function(predef_sl, plot_path, experiment=NULL){
  cont_tab <- matrix(nrow=2, ncol=max(predef_sl$df.groups), data = 0)
  t0 <- predef_sl[predef_sl$cell_id$split==0, ]
  tT <- predef_sl[predef_sl$cell_id$split==max(predef_sl$cell_id$split), ]
  for (i in 1:nrow(t0)){
    cont_tab[1,t0[i,]$df.groups] = cont_tab[1,t0[i,]$df.groups] + 1
  }
  for (i in 1:nrow(tT)){
    cont_tab[2,tT[i,]$df.groups] = cont_tab[2,tT[i,]$df.groups] + 1
  }
  
  chi2 <- chisq.test(cont_tab)
  capture.output(chi2, file = paste(plot_path, experiment, 'chi2.txt', sep = ''))
  return(chi2$p.value)
}

predef_state_transitions <- function(predef_sl, seperated_steps=F){
  # Generates an M x M state transition matrix from recorded state locations in `predef_sl`
  #
  # Parameters
  # ----------
  # predef_sl : data.frame
  #   N*T x 4 array, plate, Well.Xy, cell_id$cell, cell_id$split
  #   where N is the number of cells, T is the number of timesteps
  # separated_steps : boolean.
  #   if TRUE, returns list of `predef_stm`, indexed by timestep, each only listing transitions
  #   from a single timestep
  #
  # Returns
  # -------
  # predef_stm : data.frame.
  #   M x M transition matrix of transitions between predefined states, as provided in predef_sl
  #   Sums to N * T
  
  l <- list()
  
  for (t in min(predef_sl$cell_id$split)+1:max(predef_sl$cell_id$split)){
    state_trans_mat <- data.frame(matrix(data = 0, ncol = max(predef_sl$df.groups), nrow = max(predef_sl$df.groups)))
    
    t0_state <- predef_sl[predef_sl$cell_id$split == t-1,]
    t1_state <- predef_sl[predef_sl$cell_id$split == t,]
    for (i in 1:nrow(t0_state)){  
      cell <- t0_state[i,]
      initial <- t0_state[i,]$df.groups
      future_cell <- t1_state[ (t1_state$cell_id$cell == cell$cell_id$cell & t1_state$plate == cell$plate & t1_state$Well.XY == cell$Well.XY), ]
      final <- future_cell$df.groups
      
      state_trans_mat[initial, final] <- state_trans_mat[initial, final] + 1
      
    }
    colnames(state_trans_mat) <- 1:max(predef_sl$df.groups)
    l[[t]] <- state_trans_mat
  }
  
  
  if (seperated_steps){return(l)}
  else{
    predef_stm <- data.frame(matrix(data = 0, ncol = max(predef_sl$df.groups), nrow = max(predef_sl$df.groups)))
    for (i in 1:length(l)){
      predef_stm <- predef_stm + l[[i]]
    }
    colnames(predef_stm) <- 1:max(predef_sl$df.groups)
    return(predef_stm)
  }
  
}

predef_trans_rate_mat <- function(predef_stm){
  # Generates a transition rate matrix from a matrix of state transition observations
  # p_i,j = n_ij / sum_1^K(n_ij)
  # where p_i,j is prob of i --> j, n_i,j is num oberserved i --> j
  # K is number of states == (i || j)
  #
  # Parameters
  # ----------
  # predef_stm : N x N matrix of observed transitions i --> j
  # 
  # Returns
  # -------
  # predef_trm : N x N matrix of transition rates i --> j
  #
  
  row_sum_mat = replicate(4, apply(predef_stm, 1, sum))
  predef_trm <- predef_stm / row_sum_mat
  
  return(predef_trm)
}

plot_predef_transmat <- function(df, class_num, plot_path, method = 'ward.D2', prefix = NULL){
  
  predef_sl <- predef_state_locations(df, class_num, method=method)
  sl_contingency(predef_sl, plot_path, experiment = prefix)
  predef_stm <- predef_state_transitions(predef_sl)
  predef_trm <- predef_trans_rate_mat(predef_stm)
  write.csv(predef_trm, file = paste(plot_path, prefix, 'predef_trans_rate_matrix.txt', sep=''))
  db_broken <- calc_transmat_symmetry(predef_stm, plot_path, prefix = prefix)
  
  
  require(reshape2)
  require(ggplot2)
  require(RColorBrewer)
  tile_theme <- theme(#legend.position = "none",
    #axis.text.x=element_blank(),
    #axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    legend.position="right",
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank())
  t.m <- melt(as.matrix(predef_stm)) # Var1 = t0, Var2 = t1
  f <- numeric(nrow(t.m))
  for (i in 1:nrow(t.m)){if (t.m[i,]$Var1 == t.m[i,]$Var2){f[i] = 1}}
  frames <- t.m[f == 1, c("Var1", "Var2")]
  
  p <- ggplot(t.m, aes(x = Var2, y = Var1)) + geom_tile(aes(fill = value)) + scale_fill_distiller(palette = "RdBu")
  p <- p + tile_theme + labs(x = 't1 State', y = 't0 State') + coord_fixed() # forces square
  p <- p + geom_rect(data=frames, size=1, fill=NA, colour="black", aes(xmin=Var1 - 0.5, xmax=Var1 + 0.5, ymin=Var2 - 0.5, ymax=Var2 + 0.5))
  
  if (class(db_broken) == "data.frame"){
    rframes <- data.frame(cbind(db_broken$t0_state, db_broken$t1_state))
    colnames(rframes) <- c("Var1", "Var2")
    r2frames <- data.frame(cbind(rframes$Var2, rframes$Var1))
    colnames(r2frames) <- c("Var1", "Var2")
    rframes <- rbind(rframes, r2frames)
    p <- p + geom_rect(data=rframes, size=1, fill=NA, colour="red", aes(xmin=Var2 - 0.5, xmax=Var2 + 0.5, ymin=Var1 - 0.5, ymax=Var1 + 0.5))
  }
  
  ggsave(p, file = paste(prefix, "PFA_predef_transmat", ".png", sep=''), path = plot_path, width = 4, height = 4)
  
}

holm_bonferroni <- function(pval_mat, a = 0.05){
  # make NA holm-bonferroni matrix to hold corrected pvals
  hb_corr <- matrix(NA, nrow=nrow(pval_mat), ncol=ncol(pval_mat)) 
  # count hypothesis to be tested, assumes non p-vals are NA
  p_df <- data.frame(matrix(nrow = 0, ncol = 3))
  
  for (i in 1:nrow(pval_mat)){
    for (j in 1:ncol(pval_mat)){
      # check to ensure pval is not an NA value
      if (is.finite(pval_mat[i,j]) && pval_mat[i,j] > 0){ 
        p <- pval_mat[i,j]
        p_df <- rbind(p_df, c(p, i, j))
      }
    }
  }
  colnames(p_df) <- c("p", "i", "j")
  p_df <- p_df[order(p_df$p),]
  
  
  m <- nrow(p_df)
  k <- 1
  for (val in 1:nrow(p_df)){
    p <- p_df[val,]$p
    i <- p_df[val,]$i
    j <- p_df[val,]$j
    
    p_hb <- p * (m + 1 - k)
    # if (p_hb > a){
    #   #print('p > a')
    #   return(hb_corr) 
    #   } 
    #else {
    hb_corr[i,j] <- p_hb
    k <- k + 1
   # }
  }
  
  return(hb_corr)
}

calc_transmat_symmetry <- function(transmat, plot_path, prefix = NULL){
  upper <- transmat
  upper[lower.tri(upper, diag = T)] <- 0
  upper.t <- t(upper)
  lower <- transmat
  lower[upper.tri(lower, diag = T)] <- 0
  
  diff = upper.t - lower
  
  
  p_upper = sum(upper)/(sum(upper)+sum(lower))
  
  binom_stat <- binom.test(sum(upper), sum(upper + lower), p = 0.5, alternative = "two.sided")
  capture.output(binom_stat, file = paste(plot_path, prefix, "binom_stat.txt", sep = ''))
  capture.output(p_upper, file = paste(plot_path, prefix, "prop_upper.txt", sep = ''))
  
  
  # Test pairwise states with binom test
  
  pairwise <- matrix(NA, nrow= nrow(transmat), ncol = ncol(transmat))
  count <- 0
  for (i in 1:nrow(pairwise)){
    for (j in 1:ncol(pairwise))
      if (upper.t[i,j] > 0){
        b <- binom.test(upper.t[i,j], (upper.t[i,j]+lower[i,j]))
        pairwise[i,j] <- b$p.value
        count <- count + 1
      }
  }
  pairwise_hb <- holm_bonferroni(pairwise, a = 0.05) # HB correction, alpha = 0.05
  db_broken <- data.frame(matrix(NA, nrow=0, ncol = 3))
  for (i in 1:nrow(pairwise_hb)){
    for (j in 1:ncol(pairwise_hb)){
      if (is.finite(pairwise_hb[i,j]) && pairwise_hb[i,j] < 0.05){
        db_broken <- rbind(db_broken, c(i, j, pairwise_hb[i,j]))
      }
    }
  }
  colnames(db_broken) <- c("t0_state", "t1_state", "p-value")
  
  if (nrow(db_broken) == 0){
    db_broken <- "No pairwise transitions break detailed balance"
  }
  
  write.csv(db_broken, file = paste(plot_path, prefix, "pairwise_db_broken.txt", sep = ''))
  write.csv(pairwise_hb, file = paste(plot_path, prefix, "pairwise_p_vals.txt", sep = ''))
  return(db_broken)
}


## MuSC

experiment1 = 'musc/20160626'
experiment2 = 'musc/20160701_0'
experiment3 = 'musc/20160701_1'

predef_musc <- function(experiment1, experiment2, experiment3){
  
  plot_path = "data/musc/"
  
  df_20f1 = read.csv(paste("data/", experiment1, "/fgf2_exp_motility_statistics_split_20.csv", sep = ''))
  df_20n1 = read.csv(paste("data/", experiment1, "/nofgf2_exp_motility_statistics_split_20.csv", sep = ''))
  plate1f <- rep(1, nrow(df_20f1))
  plate1n <- rep(1, nrow(df_20n1))
  df_20f2 = read.csv(paste("data/", experiment2, "/fgf2_exp_motility_statistics_split_20.csv", sep = ''))
  df_20n2 = read.csv(paste("data/", experiment2, "/nofgf2_exp_motility_statistics_split_20.csv", sep = ''))
  plate2f <- rep(2, nrow(df_20f2))
  plate2n <- rep(2, nrow(df_20n2))
  df_20f3 = read.csv(paste("data/", experiment3, "/fgf2_exp_motility_statistics_split_20.csv", sep = ''))
  df_20n3 = read.csv(paste("data/", experiment3, "/nofgf2_exp_motility_statistics_split_20.csv", sep = ''))
  plate3f <- rep(3, nrow(df_20f3))
  plate3n <- rep(3, nrow(df_20n3))
  
  comb_df20f <- rbind(df_20f1, df_20f2, df_20f3)
  plate <- c(plate1f, plate2f, plate3f)
  comb_df20f <- data.frame(plate, comb_df20f)
  comb_df20n <- rbind(df_20n1, df_20n2, df_20n3)
  plate <- c(plate1n, plate2n, plate3n)
  comb_df20n <- data.frame(plate, comb_df20n)
  
  plot_predef_transmat(comb_df20f, class_num = 4, plot_path, prefix = "MuSC_FGF2_HClust_", method='ward.D')
  plot_predef_transmat(comb_df20n, class_num = 4, plot_path, prefix = "MuSC_noFGF2_HClust_", method='ward.D')
  
}

predef_musc(experiment1, experiment2, experiment3)

## Myoblast

predef_myoblast <- function(experiment1, experiment2){
  
  plot_path = 'data/myoblast/'
  
  df_20f1 = read.csv(paste("data/", experiment1, "/fgf2_exp_motility_statistics_split_20.csv", sep = ''))
  df_20n1 = read.csv(paste("data/", experiment1, "/nofgf2_exp_motility_statistics_split_20.csv", sep = ''))
  plate1f <- rep(1, nrow(df_20f1))
  plate1n <- rep(1, nrow(df_20n1))
  df_20f2 = read.csv(paste("data/", experiment2, "/fgf2_exp_motility_statistics_split_20.csv", sep = ''))
  df_20n2 = read.csv(paste("data/", experiment2, "/nofgf2_exp_motility_statistics_split_20.csv", sep = ''))
  plate2f <- rep(2, nrow(df_20f2))
  plate2n <- rep(2, nrow(df_20n2))

  comb_df20f <- rbind(df_20f1, df_20f2)
  plate <- c(plate1f, plate2f)
  comb_df20f <- data.frame(plate, comb_df20f)
  comb_df20n <- rbind(df_20n1, df_20n2)
  plate <- c(plate1n, plate2n)
  comb_df20n <- data.frame(plate, comb_df20n)
  
  plot_predef_transmat(comb_df20f, class_num = 2, plot_path, prefix = "Myoblast_FGF2_HClust_", method='ward.D2')
  plot_predef_transmat(comb_df20n, class_num = 2, plot_path, prefix = "Myoblast_noFGF2_HClust_", method='ward.D2')
}

myo_exp1 <- "myoblast/20160623"
myo_exp2 <- "myoblast/20160720"

predef_myoblast(myo_exp1, myo_exp2)

## MEF

plot_path = 'data/mef/'
predef_mef <- function(myc_exp, myc2_exp, wt_exp, wt2_exp, wt3_exp, wt4_exp, tau = 20, plot_path){
  
  myc_df <- read.csv(paste('data/', myc_exp, '/', 'exp_motility_statistics_split_', tau, '.csv', sep =''))
  myc2_df <- read.csv(paste('data/', myc2_exp, '/', 'exp_motility_statistics_split_', tau, '.csv', sep =''))
  wt_df <- read.csv(paste('data/', wt_exp, '/', 'exp_motility_statistics_split_', tau, '.csv', sep =''))
  wt2_df <- read.csv(paste('data/', wt2_exp, '/', 'exp_motility_statistics_split_', tau, '.csv', sep =''))
  wt3_df <- read.csv(paste('data/', wt3_exp, '/', 'exp_motility_statistics_split_', tau, '.csv', sep =''))
  wt4_df <- read.csv(paste('data/', wt4_exp, '/', 'exp_motility_statistics_split_', tau, '.csv', sep =''))
  
  
  comb_myc_df <- rbind(myc_df, myc2_df)
  plate <- c(rep(1,nrow(myc_df)), rep(2, nrow(myc2_df)))
  comb_myc_df <- data.frame(plate, comb_myc_df)
  
  comb_wt_df <- rbind(wt_df, wt2_df, wt3_df, wt4_df)
  plate <- c( rep(1, nrow(wt_df)), rep(2, nrow(wt2_df)), rep(3, nrow(wt3_df)), rep(4, nrow(wt4_df)) )
  comb_wt_df <- data.frame(plate, comb_wt_df)
  
  plot_predef_transmat(comb_myc_df, class_num = 3, plot_path, prefix = "MEF_MycRas_HClust_", method='ward.D2')
  plot_predef_transmat(comb_wt_df, class_num = 3, plot_path, prefix = "MEF_WT_HClust_", method='ward.D2')
}

myc_exp = "mef/mycras/20160917"
myc2_exp = "mef/mycras/20160918"
wt_exp = "mef/wt/20160711"
wt2_exp = "mef/wt/20160925_0"
wt3_exp = "mef/wt/20160925_1"
wt4_exp = "mef/wt/20160927"

predef_mef(myc_exp, myc2_exp, wt_exp, wt2_exp, wt3_exp, wt4_exp, tau = 20, plot_path)

## Simulations

pwr2pwr <- read.csv('data/sims/pwr2pwr_split_20.csv')
plate <- rep(1, nrow(pwr2pwr))
pwr2pwr <- data.frame(plate, pwr2pwr)
rw2rw <- read.csv('data/sims/rw2rw_split_20.csv')
plate <- rep(1, nrow(rw2rw))
rw2rw <- data.frame(plate, rw2rw)
power2rw <- read.csv('data/sims/pwr2urw_split_20.csv')
plate <- rep(1, nrow(power2rw))
power2rw <- data.frame(plate, power2rw)



plot_path = "data/sims/"
plot_predef_transmat(pwr2pwr, class_num = 3, plot_path, prefix = "pwr2pwr_HClust_")
plot_predef_transmat(rw2rw, class_num = 4, plot_path, prefix = "rw2rw_HClust_")
plot_predef_transmat(power2rw, class_num = 4, plot_path, prefix = "pwr2rw_HClust_")


## Consider all transition combinations 

predef_state_transitions_all <- function(predef_sl){
  
  state_trans_mat <- data.frame(matrix(data = 0, ncol = max(predef_sl$df.groups), nrow = max(predef_sl$df.groups)))
  
  for (j in 0:(max(predef_sl$cell_id$split)-1)){
    
    t0_state <- predef_sl[predef_sl$cell_id$split == j,]
    t1_state <- predef_sl[predef_sl$cell_id$split == j+1,]
    
    for (i in 1:nrow(t0_state)){
      cell <- t0_state[i,]
      initial <- t0_state[i,]$df.groups
      future_cell <- t1_state[ (t1_state$cell_id$cell == cell$cell_id$cell & t1_state$plate == cell$plate & t1_state$Well.XY == cell$Well.XY), ]
      
      final <- future_cell$df.groups
      state_trans_mat[initial, final] <- state_trans_mat[initial, final] + 1
      
    }
  }
  colnames(state_trans_mat) <- 1:max(predef_sl$df.groups)
  return(state_trans_mat)
  
}

plot_predef_transmat_all <- function(df, class_num, plot_path, prefix = NULL, method = 'ward.D2'){
  
  predef_sl <- predef_state_locations(df, class_num, method = method)
  predef_stm <- predef_state_transitions_all(predef_sl)
  db_broken <- calc_transmat_symmetry(predef_stm, plot_path, prefix = prefix)
  
  
  require(reshape2)
  require(ggplot2)
  require(RColorBrewer)
  tile_theme <- theme(#legend.position = "none",
    #axis.text.x=element_blank(),
    #axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    legend.position="right",
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank())
  predef_stml <- log(predef_stm)
  t.m <- melt(as.matrix(predef_stml)) # Var1 = t0, Var2 = t1
  f <- numeric(nrow(t.m))
  for (i in 1:nrow(t.m)){if (t.m[i,]$Var1 == t.m[i,]$Var2){f[i] = 1}}
  frames <- t.m[f == 1, c("Var1", "Var2")]
  p <- ggplot(t.m, aes(x = Var2, y = Var1)) + geom_tile(aes(fill = value)) + scale_fill_distiller(palette = "Spectral")
  p <- p + tile_theme + labs(x = 't1 State', y = 't0 State') + coord_fixed() # forces square
  p <- p + geom_rect(data=frames, size=1, fill=NA, colour="black", aes(xmin=Var1 - 0.5, xmax=Var1 + 0.5, ymin=Var2 - 0.5, ymax=Var2 + 0.5))
  
  if (class(db_broken) == "data.frame"){
    rframes <- data.frame(cbind(db_broken$t0_state, db_broken$t1_state))
    colnames(rframes) <- c("Var1", "Var2")
    r2frames <- data.frame(cbind(rframes$Var2, rframes$Var1))
    colnames(r2frames) <- c("Var1", "Var2")
    rframes <- rbind(rframes, r2frames)
    p <- p + geom_rect(data=rframes, size=1, fill=NA, colour="red", aes(xmin=Var2 - 0.5, xmax=Var2 + 0.5, ymin=Var1 - 0.5, ymax=Var1 + 0.5))
  }
  
  
  ggsave(p, file = paste(prefix, "PFA_predef_transmat", ".png", sep=''), path = plot_path, width = 4, height = 4)
  
}



## MuSC
plot_path = "data/musc/"
plot_predef_transmat_all(comb_df20f, class_num = 4, plot_path, prefix = "MuSC_FGF2_HClust_", method = 'ward.D')
plot_predef_transmat_all(comb_df20n, class_num = 4, plot_path, prefix = "MuSC_noFGF2_HClust_", method = 'ward.D')

## Myoblast
plot_path = "data/myoblast/"
plot_predef_transmat_all(comb_myo20f, class_num = 2, plot_path, prefix = "Myoblast_FGF2_HClust_", method = 'ward.D2')
plot_predef_transmat_all(comb_myo20n, class_num = 2, plot_path, prefix = "Myoblast_noFGF2_HClust_", method = 'ward.D2')

## MEF
plot_path = "data/mef/"
plot_predef_transmat_all(comb_myc_df, class_num = 3, plot_path, prefix = "MEF_MycRas_HClust_")
plot_predef_transmat_all(comb_wt_df, class_num = 3, plot_path, prefix = "MEF_WT_HClust_")

## Simulations

pwr2rw2fbm <- read.csv('data/split/sims/pwr2rw2fbm_ext_stat.csv')
plate <- rep(1, nrow(pwr2rw2fbm))
pwr2rw2fbm <- data.frame(plate, pwr2rw2fbm)

plot_path = "data/split/sims/"
plot_predef_transmat_all(pwr2rw2fbm, class_num = 3, plot_path, prefix = "pwr2rw2fbm_HClust_")


