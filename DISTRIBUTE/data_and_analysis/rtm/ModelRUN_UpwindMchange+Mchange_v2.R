#--------------------------------------------------------------------#
#              Vapor transport model -- model run (DIY)              #
#             -----------------------------------------              #
# T. Kukla (Stanford Univ. 2018)                                     #
#--------------------------------------------------------------------#
library(rootSolve)
library(pracma)
library(data.table)
library(MASS)
library(GGally)

rm(list=ls())  # clear global environment

# This is the main script to run the model with user-prescribed inputs
# Decisions required to initialize the simulation should be defined here
# under "CONSTRUCT DOMAIN" and "CLIMATE CONDITIONS"
# BEFORE RUNNING make sure you change the working directory 
# AND load required packages outlined in the README file
# ---- Running the entire script through will generate example model output
#      and plot the output in 4 panels in the plots viewer in Rstudio

#... identify working directory (must have other model scripts)
setwd('/DISTRIBUTE/data_and_analysis/rtm')   # [CHANGE!] to correct dir

#... bring in the function that solves the model
source('ModelSolve_fxn.R')

#... FIRST a potentially useful function:
# if you need to switch from the more common [mm yr-1] units to required
# [kg m-2 s-1] units for potential evapotranspiration, this is the 
# conversion function 
flux_conv <- function(flux){
  secs_per_min <- 60
  min_per_hr <- 60
  hr_per_day <- 24
  day_per_yr <- 365
  
  secs_per_yr <- secs_per_min * min_per_hr * hr_per_day * day_per_yr
  
  flux * (1/secs_per_yr)
}


## 1: CONSTRUCT DOMAIN ---------------------------------------------- 
n <- 600            # number of nodes in the domain
L <- 60e6            # [m] length of domain
kL <- L/1e3         # [km] length of domain
# source('build_topo_fxn.R')  # provides a function to build your domain--set equal to 'z'
z <- rep(0, n)      # flat domain
# set land mask based on land fraction (assumes all land is at the end of domain)
landbuffer <- 10e6    # [m] length of additional land at end of domain that we don't worry about
landClim <- 3e6       # [m] length of land we care about
lfrac <- (landClim + landbuffer) / L
lmask <- c(rep(0, (1-lfrac)*length(z)), rep(1, lfrac*length(z)))
# get lmask to length n
while(length(lmask) != n){
  if(length(lmask) < n){
    lmask[length(lmask)+1] <- lmask[length(lmask)]
  } else{
    lmask <- lmask[1:n]
  }
}
# ... where to record results 
gridstep <- L/n         # meters distance per 1-D grid cell
gridx <- seq(gridstep,L, by=gridstep)
coast_x <-  which.min(abs(gridx - (L - (landClim + landbuffer)))) + 1 # added small buffer
mid_x <- which.min(abs(gridx - (L - (landClim/2 + landbuffer))))
inland_x <- which.min(abs(gridx - (L - ( landbuffer))))
dx_sites <- c("coast" = coast_x, "middle" = mid_x, "inland" = inland_x)
# determine where you reach convergence zone (the last X fraction)
linstep = gridstep/L    # [m] km * meter conversion over L
csteps = 100
conv_zone_full <- scales::rescale(qexp(gridx/(L+gridstep)), to=c(landbuffer/L, 1))
conv_zone <- conv_zone_full[seq(1, length(conv_zone_full), length=csteps)]
# conv_zone <- seq(conv_zone_full[1], conv_zone_full[length(conv_zone_full)], 0.001)   # solution of stepping linearly
# conv_zone <- c(seq(landbuffer/L, (lfrac + landClim/L), by=linstep), pracma::logseq((lfrac + landClim/L)+linstep, 1, n = 30))



# determine how many random samples to execute
n_di <- 10


## 2: CLIMATE CONDITIONS --------------------------------------------
MAT_Kelv <- 305       # [K] mean annual temperature at 2m above sea level
PET <- flux_conv(1.5e3) # [kg m-2 s-1] Potential ET (flux_conv converts from mm/yr to kg/m2s)
# res.time <- 10*86400  # [s] for computing pw_initial (86400 is s/day)
res.time <- 8*86400  # [s] for computing pw_initial (86400 is s/day)
Peclet <- 100         # peclet number - ratio of advection to diffusion
pw_initial <- PET * res.time      # [kg m-2] initial vapor content (column integrated)
rh_0 <- 0.8           # relative humidity of moisture source
# u_vel <- 10           # [m s-1] advective velocity
u_vel <- 11           # [m s-1] advective velocity
trnsp <- 0.7         # transpired fraction of ET (see Good et al., 2015)
myDI.conv <- seq(0.25, 0.9, length=n_di)           # starting dryness index
myOMEGA <- 15        # prescribed recycling efficiency parameter
# update transpiration for land (zero over ocean)
trnsp <- ifelse(lmask == 1, trnsp, 0)
# update DI for convergence zone 
DI_background <- 1
# tune d18O_e to start in steady state
diff_from_min12 <- 8.96 # [per mille] difference from initial -12 vapor



# -- SHOULD DI and UPWIND CONV be correlated?? ------------------------
correlate_DI_upwind <- FALSE
# ... from https://www.r-bloggers.com/2021/05/how-to-generate-correlated-data-in-r/
if(correlate_DI_upwind == TRUE){
  # create the variance covariance matrix
  sigma<-rbind(c(3,-0.8), c(-0.8,1))
  # create the mean vector
  mu<-c(10, 5) 
  # generate the multivariate normal distribution
  dfx<-as.data.frame(mvrnorm(n=nsamples, mu=mu, Sigma=sigma))
  ggpairs(dfx)
  # rescale
  myDI.conv <- scales::rescale(dfx$V1, to=c(min(myDI.conv), max(myDI.conv)))
  conv_zone <- scales::rescale(dfx$V2, to=c(min(conv_zone), 0.66))
  plot(myDI.conv, conv_zone)
}


## 3: RUN MODEL -----------------------------------------------------
# run a control simulation (DI = 1 whole way through)
myDI <- rep(DI_background, n)
df.ctrl <- VTSolve()

# ... loop through conv_zone placements
idx.track <- 0
for(cz in 1:length(conv_zone)){
  for(i in 1:n_di){
    idx.track <- idx.track + 1
    all.iters <- length(conv_zone) * n_di
    print(paste("Now solving ", idx.track, " of ", all.iters, " -------------------", sep=''))
    # select conv_zone_placement
    conv_zone_frac <- conv_zone[cz] # sample(conv_zone, size=1)
    thisDI.conv <- myDI.conv[i]  # sample(myDI.conv, size=1)
    M_factorChange <- DI_background/thisDI.conv
    # make DI vector
    myDI <- c(rep(DI_background, (1-conv_zone_frac)*length(z)), rep(thisDI.conv, conv_zone_frac*length(z)))
    # get DI to length n
    while(length(myDI) != n){
      if(length(myDI) < n){
        myDI[length(myDI)+1] <- myDI[length(myDI)]
      } else{
        myDI <- myDI[1:n]
      }
    }
    
    # ... RUN
    df <- VTSolve() 
    
    # get the change in d18O
    for(ldx in 1:length(dx_sites)){
      this_dx <- dx_sites[ldx]
      this_dist <- gridx[this_dx]
      upwind_dist_km <- (L * conv_zone_frac/1e3) - ((L-this_dist)/1e3) # distance upwind of the specific site
      dist_inland_km <- (gridx[this_dx] - (L-L*lfrac))/1e3
      site_Dd18O <- df$d18_P[this_dx] - df.ctrl$d18_P[this_dx] 
      site_DW <- df$W[this_dx] - df.ctrl$W[this_dx]
      site_DP <- df$P[this_dx] - df.ctrl$P[this_dx]
      if((idx.track==1) & (ldx==1)){
        outdf <- as.data.table(cbind(upwind_dist_km, site_Dd18O, site_DW, site_DP, dist_inland_km, M_factorChange))
        colnames(outdf) <- c('upwind_dist_km', 'site_Dd18O', 'site_DW', 'site_DP', 'dist_inland_km', 'M_factorChange')
        outdf$site_region <- names(dx_sites[ldx])
      } else{
        tdf <- as.data.table(cbind(upwind_dist_km, site_Dd18O, site_DW, site_DP, dist_inland_km, M_factorChange))
        colnames(tdf) <- c('upwind_dist_km', 'site_Dd18O', 'site_DW', 'site_DP', 'dist_inland_km', 'M_factorChange')
        tdf$site_region <- names(dx_sites[ldx])
        outdf <- as.data.table(rbind(outdf, tdf))
      }
    }
  }
}





# ------------------------------------------
fn <- "upwindMchange+MfactorChange_results_randSamples_largerDd18.RDS"
fncsv <- "upwindMchange+MfactorChange_randSamples_results_largerDd18.csv"
# saveRDS(outdf, fn)
# write.csv(outdf, fncsv)
# *************************************************************************************************** #
# --------------------------------------------------------------------------------------------------- #
# *************************************************************************************************** #

