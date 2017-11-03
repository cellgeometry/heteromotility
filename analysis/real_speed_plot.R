# Plot real speeds for all cell types
require(reshape2)
require(ggplot2)

mef <- read.csv('data/mef/MEF_speeds_real_units.csv')
mef$celltype <- 'MEF'
mef$typeclass <- paste(mef$celltype, mef$Class, sep=' ')
myoblast <- read.csv('data/myoblast/Myoblast_speeds_real_units.csv')
myoblast$celltype <- 'Myoblast'
myoblast$typeclass <- paste(myoblast$celltype, myoblast$Class, sep=' ')
musc <- read.csv('data/musc/MuSC_speeds_real_units.csv')
musc$celltype <- 'MuSC'
musc$typeclass <- paste(musc$celltype, musc$Class, sep=' ')

df <- rbind(mef, myoblast, musc)
df.m <- melt(df, id.vars=c('typeclass', 'Class', 'celltype'))

# plot

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

p <- ggplot(data=df, aes(x=typeclass, y=SpeedMean, fill=typeclass)) + geom_bar(stat='identity')
p <- p + geom_errorbar(aes(ymin=SpeedMean-SpeedSE, ymax=SpeedMean+SpeedSE)) 
p <- p + labs(x='', y='Velocity [um/min]', title='Cell Type Velocities') + guides(fill=guide_legend(title="Cell Type"))
p <- p + scale_y_continuous(expand = c(0,0), limits = c(0, 1.10*max(df$SpeedMean+df$SpeedSE)))
p <- p + theme(axis.text.x=element_text(angle=45,hjust=1)) + cluster_theme

ggsave(p, path = 'data/comparisons/', filename='celltype_real_speeds.png', height = 4, width = 4)
