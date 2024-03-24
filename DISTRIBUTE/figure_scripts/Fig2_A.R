# ----------------------------------------------
#
# plot the distribution of d18O responses to 
# precession (requires output of:
# "SISAL_smoothRange_perOrbCycle.R" script)
# 
# ---
# T Kukla (UW Atmospheric Sciences, 2023)
#
# ----------------------------------------------
library(data.table)
library(ggplot2) 
library(ggthemes)
library(palinsol)
library(maps)

rm(list=ls())


# --- READ IN DATA
mainPath <- '/DISTRIBUTE/data_and_analysis/sisal_v2_d18Oranges'  # [CHANGE!] directory to shown folder ("sisla_v2_d18Oranges")
setwd(mainPath)

# ... read in dat
df <- readRDS('SISALv2_2kyrSmoothRangesINSOL_perOrbCycle.RDS')

# WHERE TO SAVE
save.here <- mainPath

# --- PLOT
# ... median line
med <- median(df[latitude > -37 & latitude < 37]$xrange)
# ... 90th percentile
p90 <- quantile(df[latitude > -37 & latitude < 37]$xrange, probs=0.9)


# [1] HISTOGRAM
p.hist <- ggplot(df[latitude > -37 & latitude < 37]) +
  coord_cartesian(ylim=c(0,0.5)) +
  # median line
  annotate(geom='segment', x=med, xend=med, y=0, yend=0.5,
           size=1, color='darkred', linetype='dashed') +
  # 90th percentile
  annotate(geom='segment', x=p90, xend=p90, y=0, yend=0.5,
           size=1, color='pink', linetype='dashed') +
  # hist and density
  geom_histogram(aes(x=xrange, y=..density..), bins=20, color='darkgray') +
  geom_density(aes(x=xrange), size=2) + 
  scale_y_continuous(name="density", expand=c(0,0)) +
  scale_x_continuous(name=expression("smoothed δ"^"18"*"O"~"range (‰)"),
                     limits=c(0,10), breaks=c(0:10)) +
  # lims(x=c(0, 3)) + 
  theme_few() +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.text.x = element_text(size=20), axis.title = element_text(size=25))

p.hist

save.fn <- "orbCycle_Hist.png"
save.name <- paste(save.here, save.fn, sep='/')
# ggsave(save.name, p.hist, width=18, height=15, units='cm')
