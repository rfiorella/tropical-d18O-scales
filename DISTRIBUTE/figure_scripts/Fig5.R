# ------------------------------------------
# 
# upwind M change model PLOT 
#
# ------------------------------------------
library(ggplot2)
library(ggthemes)
library(data.table)
library(cmocean)
library(ggh4x)
library(akima)

rm(list=ls())

setwd('/Users/tylerkukla/Documents/academic_research/myProjects/spatialScales_tropicalIsotopes/DISTRIBUTE/data_and_analysis/rtm')   # [CHANGE!] to correct path

# ... read in data
fn <- "upwindMchange+MfactorChange_results_randSamples_largerDd18.RDS"
df <- readRDS(fn)

# ... parameters copied over from run script
L <- 60e6            # [m] length of domain
# set land mask based on land fraction (assumes all land is at the end of domain)
landbuffer <- 10e6    # [m] length of additional land at end of domain that we don't worry about

# ... alphabetize the site regions
df$site_region_alpha <- "X"
df[site_region=='coast']$site_region_alpha <- "A_coast"
df[site_region=='middle']$site_region_alpha <- "B_middle"
df[site_region=='inland']$site_region_alpha <- "C_inland"

# plot decisions
colorbar = "amp"
tick_fsize = 20
lab_fsize = 24
leg_fsize=20


# ... INTERPOLATE 
n_pts <- 500
dfi <- with(df[site_region=='inland'], akima::interp(x = upwind_dist_km/1e4, y = M_factorChange, z = site_Dd18O,
                       linear = TRUE, extrap = FALSE, duplicate="mean",
                       xo = seq(1,6e4, length=n_pts)/1e4,
                       yo = seq(0,6, length=n_pts)))
dfi.p <- as.data.frame(interp2xyz(dfi))
dfi.p$x <- dfi.p$x * 1e4


# --- PLOT
library(metR)

save.here <- ''  # where to save resulting panel

p.cont <- ggplot(dfi.p) + 
  coord_cartesian(ylim=c(1.15,3.99), xlim=c(50,5e4)) +
  geom_contour_fill(aes(x=x, y=y, z=abs(z)), breaks=c(0:6)) +
  scale_x_log10("distance to rainout change (km)", expand=c(0,0)) + 
  scale_y_continuous("rainout intensity factor increase", expand=c(0,0),
                     breaks=c(1.5,2,2.5,3,3.5)) +
  scale_fill_cmocean("", name="rain", direction=-1) +
  theme_few() +
  theme(axis.text = element_text(size=20), axis.title = element_text(size=25))

p.cont

save.fn <- "RTM_spaceVint.png"
save.name <- paste(save.here, save.fn, sep='/')
# ggsave(save.name, p.cont, width=24, height=16, units='cm')




# ------------

