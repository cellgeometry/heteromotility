# variance explained

musc <- read.csv('data/musc/MuSC_pca_cum_var_explained.csv')
mef <- read.csv('data/mef/MEF_pca_cum_var_explained.csv')
myoblast <- read.csv('data/myoblast/Myoblast_pca_cum_var_explained.csv')

plt_df <- data.frame(musc, mef, myoblast, 1:nrow(musc))
colnames(plt_df) <- c('MuSC', 'MEF', 'Myoblast', 'PC')

plt_df.m <- melt(plt_df, id.vars=c('PC'))

p <- ggplot(plt_df.m, aes(PC, y=value, color=variable)) + geom_freqpoly(stat='identity', size=0.8)
p <- p + labs(title='Cumulative Explained Variance by PC', x='PC', y='Cumulative Explained Variance')
p <- p + scale_y_continuous(expand=c(0,0), limits=c(0.3,1.1)) + guides(color=guide_legend(title="Cell Type"))
ggsave(p, path='data/comparisons/', filename = 'cumulative_pc_variance.png', width=4, height=3, units='in')
