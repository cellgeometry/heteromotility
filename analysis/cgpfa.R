### Satellite Cell PFA for combined three plates


pfa_state_locations <- function(df, bin_sizes){
  require(reshape2)
  require(ash)
  # PCA to generate 2 dimensional state space
  df.vals <- data.frame( scale( cbind(df[4:22], df[29:32], df[39:42], df[49:52]) ) )
  df.pca <- prcomp(df.vals)
  df.comp <- data.frame( df.pca$x[,1:2] )
  pca_frame <- cbind(df[,1:3], df.comp)
  pca_frame <- transform(pca_frame, cell_id = colsplit(cell_id, "-", names = c("cell", "split")))
  # Find PC ranges to bin space
  PC1_range <- c( floor(range(df.comp[,1])[1]), ceiling(range(df.comp[,1])[2]) )
  PC2_range <- c( floor(range(df.comp[,2])[1]), ceiling(range(df.comp[,2])[2]) )

  require(ash)
  ab <- matrix( c(PC1_range[1], PC2_range[1], PC1_range[2], PC2_range[2]), 2, 2 )
  
  # Cols = Plate, XY, Split, Cell, BinRow(Y), BinCol(X)
  # Rows = Cells
  state_locations <- data.frame(matrix(ncol = 6, nrow = nrow(df)))
  colnames(state_locations) <- c("plate", "xy", "cell", "split", "biny", "binx")
  for (i in 1:nrow(pca_frame)){
    cell <- pca_frame[i,]
    position <- cell[,4:5]
    bins <- bin2(position, ab, nbin = bin_sizes)
    bin_yx <- which(bins$nc != 0, arr.ind = T)
    state_locations[i,] <- c(cell$plate, as.character(cell$Well.XY), cell$cell_id$cell, cell$cell_id$split, bin_yx)
  }
  
  state_locations$plate <- as.numeric(state_locations$plate)
  state_locations$cell <- as.numeric(state_locations$cell)
  state_locations$split <- as.numeric(state_locations$split)
  state_locations$binx <- as.numeric(state_locations$binx)
  state_locations$biny <- as.numeric(state_locations$biny)
  
  return(state_locations)
}

pfa_state_vectors <- function(state_locations){
  ## Count State Transitions
  #
  # Returns
  # -------
  # state_vectors : data.frame.
  # N x 5 matrix, where N is the number of unit transitions observed in the timecourse
  # Columns: biny, binx, vy, vx, cells
  # binx and biny are the coordinates for each course-grained bin
  # vx and vy are the corresponding components of the mean transition vector
  # cells is the total number of cells observed in the state over time
  bin_transitions <- matrix(ncol = 4, nrow = 0)   # cols = (biny, binx, v_y, v_x)
  for (s in 0:(max(state_locations$split)-1)){
    # Make data.frames for each state time scale
    t0_state <- subset(state_locations, state_locations$split==s)
    t1_state <- subset(state_locations, state_locations$split==s+1)
    # Check to ensure all cells are in both state periods
    occupied_bins <- unique(t0_state[,c("biny","binx")])
    
    for (i in 1:nrow(occupied_bins)){
      cells <- subset(state_locations, (state_locations$biny == occupied_bins[i,1] & state_locations$binx == occupied_bins[i,2]))
      for (j in 1:nrow(cells)){
        position <- c(cells[j,]$biny, cells[j,]$binx) # position at t0
        trans <- subset(t1_state, (t1_state$plate == cells[j,1] & t1_state$xy == cells[j,2] & t1_state$cell == cells[j,3]))
        trans_pos <- c(trans$biny, trans$binx)
        v <- trans_pos - position
        print(i)
        print(j)
        while ( sum(sqrt(v**2)) > 0 ){
          # Performs path finding by moving in the direction of greater magnitude
          # If direction mags are ==; randomly chooses direction
          if ( abs(v[1]) > abs(v[2]) ){
            bin_transitions <- rbind(bin_transitions, c(position, sign(v[1]), 0))
            position <- position + c(sign(v[1]), 0)
            v <- v - c(sign(v[1]), 0)
          } else if ( abs(v[1]) == abs(v[2]) ){
            choice <- rbinom(n = 1, size = 1, prob = 0.5)
            if (choice == 0){
              bin_transitions <- rbind(bin_transitions, c(position, sign(v[1]), 0))
              position <- position + c(sign(v[1]), 0)
              v <- v - c(sign(v[1]), 0)
            } else {
              bin_transitions <- rbind(bin_transitions, c(position, 0, sign(v[2])))
              position <- position + c(0, sign(v[2]))
              v <- v - c(0, sign(v[2]))
            }
          }
          else {
            bin_transitions <- rbind(bin_transitions, c(position, 0, sign(v[2])))
            position <- position + c(0, sign(v[2]))
            v <- v - c(0, sign(v[2]))
          }
        }
      }
    }
  }
  bin_transitions <- data.frame(bin_transitions)
  colnames(bin_transitions) <- c("biny","binx","vector_y", "vector_x")
  bin_transitions$binx <- as.numeric(bin_transitions$binx)
  bin_transitions$biny <- as.numeric(bin_transitions$biny)
  bin_transitions$vector_y <- as.numeric(bin_transitions$vector_y)
  bin_transitions$vector_x <- as.numeric(bin_transitions$vector_x)
  
  ## Find State Vectors
  
  state_vectors <- data.frame(matrix(ncol=5, nrow=0))
  
  i = 1
  y_bins <- sort(unique(bin_transitions$biny))
  
  for (y in y_bins){
    y_df <- subset(bin_transitions, bin_transitions$biny == y)
    x_bins <- sort(unique(y_df$binx))
    for (x in x_bins){
      xy_df <- subset(y_df, y_df$binx == x)
      v_mean <- c(mean(xy_df$vector_y), mean(xy_df$vector_x))
      cell_num <- nrow(xy_df)
      state_vectors <- rbind(state_vectors, c(y, x, v_mean, cell_num))
    }
  }
  colnames(state_vectors) <- c("biny", "binx", "v_y", "v_x", "cells")
  return(state_vectors) 
}

pfa_vector_distribution <- function(state_locations){
  
  mat_rows = sum(state_locations$split == 1)+sum(state_locations$split==2)
  vector_dist <- matrix(nrow=0, ncol = 3)
  for (s in 0:(max(state_locations$split)-1)){
    t0_state <- subset(state_locations, state_locations$split==s)
    t1_state <- subset(state_locations, state_locations$split==s+1)
    
    for ( i in 1:nrow(t0_state) ){
      v <- c((t1_state[i,]$biny - t0_state[i,]$biny),(t1_state[i,]$binx - t0_state[i,]$binx))
      v_mag <- sqrt(sum(v**2))
      vector_dist <- rbind(vector_dist, c(v, v_mag))
    }
  }
  vd <- data.frame(vector_dist)
  colnames(vd) <- c("v_y", "v_x", "v_mag")
  return(vd)
}

sem <- function(x){sd(x)/sqrt(length(x))}

vector_distribution_stats <- function(state_locations, out_path, prefix = NULL){
  require(moments)
  vd <- pfa_vector_distribution(state_locations)
  write.csv(vd$v_mag, file = paste(out_path, prefix, "v_mag_vals.csv", sep=''), row.names = F)
  mean_mag <- mean(vd$v_mag)
  var_mag <- var(vd$v_mag)
  skew_mag <- skewness(vd$v_mag)
  kurt_mag <- kurtosis(vd$v_mag)
  sem_mag <- sem(vd$v_mag)
  n <- nrow(vd)
  # prob of flux for a given cell
  p_flux <- sum(vd$v_mag > 0)/nrow(vd)
  # directedness vector
  dir_v <- c(mean(vd$v_x), mean(vd$v_y))
  mag_dir <- sqrt(sum(dir_v**2))
  r <- data.frame(mean_mag, var_mag, skew_mag, kurt_mag, sem_mag, n, p_flux, dir_v[1], dir_v[2], mag_dir)
  colnames(r) <- c('mean_mag', 'var_mag', 'skew_mag', 'kurt_mag', 'sem_mag', 'n', 'p_flux', 'dir_vx', 'dir_vy', 'mag_dir')
  write.csv(r, file = paste(out_path, prefix, "vector_dist_stats.csv", sep =''), row.names = F)
  
  p <- ggplot(vd, aes(v_mag)) + geom_density()
  ggsave(p, filename = paste(prefix, "vector_dist_density.png", sep=''), path = out_path, width = 3, height = 3)
}

output_state_vectors <- function(state_vectors, output_file){
  output <- data.frame(cbind(state_vectors$binx, state_vectors$biny, state_vectors$v_x, state_vectors$v_y, state_vectors$cells))
  colnames(output) <- c("x", "y", "v_x", "v_y", "cells")
  write.csv(output, file = output_file, row.names = F)
}

save_state_vectors <- function(df, bin_sizes, output_file){
  state_locations <- pfa_state_locations(df, bin_sizes)
  state_vectors <- pfa_state_vectors(state_locations)
  output_state_vectors(state_vectors, output_file)
}

pfa_flux_plot <- function(df, bin_sizes, output_path, experiment = NULL){
  require(ggplot2)
  
  state_locations <- pfa_state_locations(df, bin_sizes)
  state_vectors <- pfa_state_vectors(state_locations)
  
  ## Vector field plot
  v_plot <- ggplot(state_vectors, aes(x = binx, y = biny)) + geom_segment(aes(xend = binx+v_x, yend = biny+v_y), arrow = arrow(length = unit(0.1, "cm")))
  v_plot <- v_plot + labs(title = "Probability Flux (tau0 - tau1)", x = "PC1", y = "PC2")
  ggsave(v_plot, file = paste(experiment, "pfa_flux_vectors.png", sep = ''), path = output_path, width = 4, height = 4)
  output_file = paste(output_path, experiment, "state_vectors.csv", sep = '')
  output_state_vectors(state_vectors, output_file)
  vector_distribution_stats(state_locations, output_path, prefix = experiment)
  return(v_plot)
}

## MuSCs

experiment1 = 'musc/20160626'
experiment2 = 'musc/20160701_0'
experiment3 = 'musc/20160701_1'

pfa_musc <- function(experiment1, experiment2, experiment3, plot_path='data/musc/', bin_sizes=c(15,15), tau=20){
  
  df_20f1 = read.csv(paste("data/", experiment1, "/fgf2_exp_motility_statistics_split_", tau, ".csv", sep = ''))
  df_20n1 = read.csv(paste("data/", experiment1, "/nofgf2_exp_motility_statistics_split_", tau, ".csv", sep = ''))
  plate1f <- rep(1, nrow(df_20f1))
  plate1n <- rep(1, nrow(df_20n1))
  df_20f2 = read.csv(paste("data/", experiment2, "/fgf2_exp_motility_statistics_split_", tau, ".csv", sep = ''))
  df_20n2 = read.csv(paste("data/", experiment2, "/nofgf2_exp_motility_statistics_split_", tau, ".csv", sep = ''))
  plate2f <- rep(2, nrow(df_20f2))
  plate2n <- rep(2, nrow(df_20n2))
  df_20f3 = read.csv(paste("data/", experiment3, "/fgf2_exp_motility_statistics_split_", tau, ".csv", sep = ''))
  df_20n3 = read.csv(paste("data/", experiment3, "/nofgf2_exp_motility_statistics_split_", tau, ".csv", sep = ''))
  plate3f <- rep(3, nrow(df_20f3))
  plate3n <- rep(3, nrow(df_20n3))
  
  comb_df20f <- rbind(df_20f1, df_20f2, df_20f3)
  plate <- c(plate1f, plate2f, plate3f)
  comb_df20f <- data.frame(plate, comb_df20f)
  comb_df20n <- rbind(df_20n1, df_20n2, df_20n3)
  plate <- c(plate1n, plate2n, plate3n)
  comb_df20n <- data.frame(plate, comb_df20n)
  
  v_comb_df20f <- pfa_flux_plot(comb_df20f, bin_sizes = bin_sizes, output_path = plot_path, experiment = paste("MuSC_fgf2_", tau, '_', 'b', bin_sizes[1], '_', sep=''))
  v_comb_df20n <- pfa_flux_plot(comb_df20n, bin_sizes = bin_sizes, output_path = plot_path, experiment = paste("MuSC_nofgf2_", tau, '_', 'b', bin_sizes[1], '_', sep=''))
  
}

bin_sizes_all = rbind(c(3,3), c(5,5), c(10,10), c(15,15), c(20,20), c(30,30))
for (tau in c(20, 25, 30)){
  for (j in 1:nrow(bin_sizes_all)){
    pfa_musc(experiment1, experiment2, experiment3, tau=tau, bin_sizes=bin_sizes_all[j,])
  }
}

## Myoblasts


pfa_myoblast <- function(myo_exp1, myo_exp2, tau = 20, bin_sizes=c(15,15), plot_path){
  myo1_f <- read.csv(paste('data/', myo_exp1, '/', 'fgf2_exp_motility_statistics_split_', tau, '.csv', sep = ''))
  myo1_n <- read.csv(paste('data/', myo_exp1, '/', 'nofgf2_exp_motility_statistics_split_', tau, '.csv', sep = ''))
  myo2_f <- read.csv(paste('data/', myo_exp2, '/', 'fgf2_exp_motility_statistics_split_', tau, '.csv', sep = ''))
  myo2_n <- read.csv(paste('data/', myo_exp2, '/', 'nofgf2_exp_motility_statistics_split_', tau, '.csv', sep = ''))
  
  comb_myo20f <- rbind(myo1_f, myo2_f)
  plate <- c(rep(1, nrow(myo1_f)), rep(2, nrow(myo2_f)))
  comb_myo20f <- data.frame(plate, comb_myo20f)
  
  comb_myo20n <- rbind(myo1_n, myo2_n)
  plate <- c(rep(1, nrow(myo1_n)), rep(2, nrow(myo2_n)))
  comb_myo20n <- data.frame(plate, comb_myo20n)
  
  v_myo20f <- pfa_flux_plot(comb_myo20f, bin_sizes = bin_sizes, plot_path, experiment = paste("Myoblast_FGF2_", tau, '_b', bin_sizes[1], '_', sep=''))
  v_myo20n <- pfa_flux_plot(comb_myo20n, bin_sizes = bin_sizes, plot_path, experiment = paste("Myoblast_noFGF2_", tau, '_b', bin_sizes[1], '_', sep=''))
  
}

myo_exp1 <- "myoblast/20160623"
myo_exp2 <- "myoblast/20160720"
plot_path = "data/myoblast/"

bin_sizes_all = rbind(c(3,3), c(5,5), c(10,10), c(15,15), c(20,20), c(30,30))
for (tau in c(20, 25, 30)){
  for (j in 1:nrow(bin_sizes_all)){
    pfa_myoblast(myo_exp1, myo_exp2, tau=tau, bin_sizes = bin_sizes_all[j,], plot_path=plot_path)
  }
}

## MEFs

plot_path = 'data/mef/'
pfa_mycwt_combined <- function(myc_exp, myc2_exp, wt_exp, wt2_exp, wt3_exp, wt4_exp, tau = 20, bin_sizes=c(15,15), plot_path){
  
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
  
  v_myc <- pfa_flux_plot(df = comb_myc_df, bin_sizes = bin_sizes, plot_path, experiment = paste("MEF_MycRas_", tau, '_b', bin_sizes[1], '_', sep=''))
  v_wt <- pfa_flux_plot(df = comb_wt_df, bin_sizes = bin_sizes, plot_path, experiment = paste("MEF_WT_", tau, '_b', bin_sizes[1], '_', sep=''))
}

myc_exp = "mef/mycras/20160917"
myc2_exp = "mef/mycras/20160918"
wt_exp = "mef/wt/20160711"
wt2_exp = "mef/wt/20160925_0"
wt3_exp = "mef/wt/20160925_1"
wt4_exp = "mef/wt/20160927"


bin_sizes_all = rbind(c(3,3), c(5,5), c(10,10), c(15,15), c(20,20), c(30,30))
for (tau in c(20, 25, 30)){
  for (j in 1:nrow(bin_sizes_all)){
    pfa_mycwt_combined(myc_exp, myc2_exp, wt_exp, wt2_exp, wt3_exp, wt4_exp, tau=tau, bin_sizes=bin_sizes_all[j,], plot_path=plot_path)
  }
}


### Combined Simulations
plot_path = "data/sims/"

fbm2rw <- read.csv('data/sims/fbm2urw_split_20.csv')
plate <- rep(1, nrow(fbm2rw))
fbm2rw <- data.frame(plate, fbm2rw)
power2fbm <- read.csv('data/sims/pwr2fbm_split_20.csv')
plate <- rep(1, nrow(power2fbm))
power2fbm <- data.frame(plate, power2fbm)
power2rw <- read.csv('data/sims/pwr2urw_split_20.csv')
plate <- rep(1, nrow(power2rw))
power2rw <- data.frame(plate, power2rw)

v_fbm2rw <- pfa_flux_plot(fbm2rw, c(15,15), plot_path, experiment = 'fbm2rw_')
v_pwr2fbm <- pfa_flux_plot(power2fbm, c(15,15), plot_path, experiment = 'pwr2fbm_')

for (j in 1:nrow(bin_sizes_all)){
  v_pwr2rw <- pfa_flux_plot(power2rw, bin_sizes_all[j,], plot_path, experiment = paste('pwr2rw_',  'b', bin_sizes_all[j,1], '_', sep=''))
}
# controls

pwr2pwr <- read.csv('data/sims/pwr2pwr_split_20.csv')
plate <- rep(1, nrow(pwr2pwr))
pwr2pwr <- data.frame(plate, pwr2pwr)
rw2rw <- read.csv('data/sims/rw2rw_split_20.csv')
plate <- rep(1, nrow(rw2rw))
rw2rw <- data.frame(plate, rw2rw)

v_pwr2pwr <- pfa_flux_plot(pwr2pwr, c(15,15), plot_path, experiment = "pwr2pwr_")


v_rw2rw <- pfa_flux_plot(rw2rw, c(15,15), plot_path, experiment = "rw2rw_")

bin_sizes_all = rbind(c(3,3), c(5,5), c(10,10), c(15,15), c(20,20), c(30,30))
for (j in 1:nrow(bin_sizes_all)){
  v_rw2rw <- pfa_flux_plot(rw2rw, bin_sizes_all[j,], plot_path, experiment = paste("rw2rw_", 'b', bin_sizes_all[j,1], '_', sep=''))
}

## Compare MEF and MuSC vector distributions

compare_vdists <- function(df1, df2, plot_path, exp1 = NULL, exp2 = NULL){
  vd1 <- pfa_vector_distribution(pfa_state_locations(df1, c(15,15)))
  vd2 <- pfa_vector_distribution(pfa_state_locations(df2, c(15,15)))
  t <- t.test(vd1$v_mag, vd2$v_mag, alternative = 'two.sided', var.equal = F, conf.level = 0.95)
  capture.output(t, file = paste(plot_path, exp1, exp2, 'vmag_ttest.txt', sep = ''))
  return(t)
}

plot_path = "data/split/vdist_comparisons/"

compare_vdists(comb_df20f, comb_myc_df, plot_path, exp1 = "MuSC_fgf2_", exp2 = "MEF_MycRas_")
compare_vdists(comb_df20f, comb_wt_df, plot_path, exp1 = "MuSC_fgf2_", exp2 = "MEF_WT_")

compare_vdists(comb_df20n, comb_myc_df, plot_path, exp1 = "MuSC_nofgf2_", exp2 = "MEF_MycRas_")
compare_vdists(comb_df20n, comb_wt_df, plot_path, exp1 = "MuSC_nofgf2_", exp2 = "MEF_WT_")

compare_vdists(comb_df20f, comb_df20n, plot_path, exp1 = "MuSC_fgf2_", exp2 = "MuSC_nofgf2_")
compare_vdists(comb_myc_df, comb_wt_df, plot_path, exp1 = "MEF_MycRas_", exp2 = "MEF_WT_")

compare_vdists(power2rw, pwr2pwr, plot_path, exp1 = "pwr2rw_", exp2 = "pwr2pwr_")
compare_vdists(power2rw, rw2rw, plot_path, exp1 = "pwr2rw_", exp2 = "rw2rw_")
compare_vdists(pwr2pwr, rw2rw, plot_path, exp1 = "pwr2pwr_", exp2 = "rw2rw_")




