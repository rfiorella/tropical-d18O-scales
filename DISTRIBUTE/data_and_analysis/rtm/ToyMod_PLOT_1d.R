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

rm(list=ls())

setwd('C:/Users/tkukl/OneDrive/Documents/MonsoonEBM+isoAttenuation/rtm')

# ... read in data
# fn <- "upwindM_change_results.RDS"
fn <- "upwindM_change_results_largerDd18.RDS"
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



# --- PLOT 
p.x <- ggplot(df) + 
  coord_cartesian(xlim =c(100, ((L-landbuffer)/1e3)-3e3), ylim=c(0,6)) + 
  annotate(geom='rect', xmin=1e3, xmax=5e3, ymin=-3, ymax=10, fill="#FBCF8E", alpha=0.4) +
  geom_line(aes(x=upwind_dist_km, y=abs(site_Dd18O), group=site_region_alpha, color=site_region_alpha), size=3) + 
  scale_color_cmocean(discrete=T, labels=c("A_coast"="coast", "B_middle" = "middle", "C_inland" = "inland"), name=colorbar,
                      start=0.2, end=1) + 
  scale_x_log10(expand=c(0,0), breaks=c(100, 1e3, 1e4 ), name="Upwind distance of rainout increase (km)",
                minor_breaks = seq(100, 1000, by = 100)) +
  annotation_logticks(sides="b") + 
  scale_y_continuous(name=expression("|δ"^"18"*"O"~"anomaly| (‰)")) +
  theme_few() +
  theme(legend.title = element_blank(), axis.text = element_text(size=tick_fsize), 
        axis.title = element_text(size=lab_fsize), legend.text = element_text(size=leg_fsize))

p.x

save.here <- 'C:/Users/tkukl/OneDrive/Documents/Talks/AGU/2023/presentation'
save.fn <- "ToyMod_3sites_largerShift.png"
save.name <- paste(save.here, save.fn, sep='/')
ggsave(save.name, p.x, width=25, height=16, units='cm')






