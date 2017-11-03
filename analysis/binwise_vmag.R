# Plot simple bar graph of binwise vector magnitudes

require(ggplot2)
require(reshape2)

# Read in transition mags (binwise)

calc_binwise_vmag <- function(tau=20, bins=15){
  
  mycras <- read.csv(paste('divergence_figs/MEF_MycRas_t', tau, '_b', bins, '_bin_vmag.csv', sep=''))
  wt <- read.csv(paste('divergence_figs/MEF_WT_t', tau, '_b', bins, '_bin_vmag.csv', sep=''))
  
  musc_fgf <- read.csv(paste('divergence_figs/MuSC_FGF2_t', tau, '_b', bins, '_bin_vmag.csv', sep=''))
  musc_nofgf <- read.csv(paste('divergence_figs/MuSC_noFGF2_t', tau, '_b', bins, '_bin_vmag.csv', sep=''))
  
  myo_fgf <- read.csv(paste('divergence_figs/Myoblast_FGF2_t', tau, '_b', bins, '_bin_vmag.csv', sep=''))
  myo_nofgf <- read.csv(paste('divergence_figs/Myoblast_noFGF2_t', tau, '_b', bins, '_bin_vmag.csv', sep=''))
  
  pwr2rw <- read.csv(paste('divergence_figs/pwr2rw_b', bins, '_bin_vmag.csv', sep = ''))
  rw2rw <- read.csv(paste('divergence_figs/rw2rw_b', bins, '_bin_vmag.csv', sep = ''))
  
  # Collect means and SEMs for plotting
  
  rep2vector <- function(vals, reps){
    # Repeats elements of vector vals by number of times in reps
    # Returns concatenated expanded vector
    newv <- numeric(0)
    for (i in 1:length(vals)){
      newv <- c(newv, rep(vals[i], reps[i]))
    }
    return(newv)
  }
  
  plt_mat <- matrix(nrow=8, ncol=2)
  
  i = 1
  for (v in list(mycras, wt, musc_fgf, musc_nofgf, myo_fgf, myo_nofgf, pwr2rw, rw2rw)){
    # sum of weighted vmags is mean v mag
    # SEM is SD of concatenated bin_vmag * cells / sqrt(sum(cells))
    sem <- sd(rep2vector(v$bin_vmag, v$cells)) / sqrt(sum(v$cells))
    plt_mat[i,] <- c(sum(v$bin_vmag_w), sem)
    i = i + 1
  }
  
  # Statistics
  
  control_t <- t.test(rep2vector(pwr2rw$bin_vmag, pwr2rw$cells), rep2vector(rw2rw$bin_vmag, rw2rw$cells))
  musc_v_wt_t <- t.test(rep2vector(musc_fgf$bin_vmag, musc_fgf$cells), rep2vector(wt$bin_vmag, wt$cells))
  musc_v_mycras_t <- t.test(rep2vector(musc_fgf$bin_vmag, musc_fgf$cells), rep2vector(mycras$bin_vmag, mycras$cells))
  
  
  
  # Plot
  
  plt_df <- data.frame(plt_mat)
  colnames(plt_df) <- c('Mean', 'SE')
  plt_df$System <- c('MEF MycRas', 'MEF WT', 'MuSC FGF2+', 'MuSC FGF2-', 'Myo. FGF2+', 'Myo. FGF2-', 'Flier to RW', 'RW only')
  plt_df$System <- factor(plt_df$System, levels=plt_df$System)
  
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
  
  p <- ggplot(data=plt_df, aes(x=System, y=Mean, fill=System)) 
  p <- p + geom_bar(stat='identity') + geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE)) 
  p <- p + labs(x='', y='Statewise Transition Mag. (cgPFA units)', title='Transition Directionality')
  p <- p + cluster_theme + scale_y_continuous(expand = c(0,0), limits = c(0, 1.10*max(plt_df$Mean+plt_df$SE)))
  p <- p + theme(axis.text.x=element_text(angle=45,hjust=1))
  
  ggsave(plot = p, filename = paste('Statewise_Transition_Directionality_t', tau, '_b', bins, '.png', sep=''), path='divergence_figs/', width = 4, height=4)
  
}
