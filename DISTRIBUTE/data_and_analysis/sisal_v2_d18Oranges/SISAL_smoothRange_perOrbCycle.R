# ----------------------------------------------
#
# loop through tropical SISALv2 data and 
# collect information about the smoothed 
# range across precession cycles
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

mainPath <- '/Users/tylerkukla/Documents/academic_research/myProjects/spatialScales_tropicalIsotopes/DISTRIBUTE/data_and_analysis/sisal_v2_d18Oranges'
setwd(mainPath)

# ************************************* #
# Define the relevant days of year for  
# computation (normal [non-leap] year)
# (after: https://cals.arizona.edu/azmet/julian.html)
Jan.start <- 1; Jan.end <- 31
Feb.start <- 32; Feb.end <- 59
Mar.start <- 60; Mar.end <- 90
Apr.start <- 91; Apr.end <- 120
May.start <- 121; May.end <- 151
Jun.start <- 152; Jun.end <- 181
Jul.start <- 182; Jul.end <- 212
Aug.start <- 213; Aug.end <- 243
Sep.start <- 244; Sep.end <- 273
Oct.start <- 274; Oct.end <- 304
Nov.start <- 305; Nov.end <- 334
Dec.start <- 335; Dec.end <- 365


# ----------------------------------------------------------------------------- 
# read in data (one at a time)
# ----------------------------------------------------------------------------- 

#... Isotopes
d13 <- as.data.table(read.csv('d13C.csv'))
d18 <- as.data.table(read.csv('d18O.csv'))

#... Sites 
site <- as.data.table(read.csv('site.csv'))

#... Entity
entity <- as.data.table(read.csv('entity.csv'))

#... Chronology
ages <- as.data.table(read.csv('original_chronology.csv'))

#... Samples
sampledf <- as.data.table(read.csv('sample.csv'))


# ----------------------------------------------------------------------------- 
# bring it all together
# ----------------------------------------------------------------------------- 

d13d18 <- merge(d13, d18, by=c('sample_id'), all=TRUE)

sample_ent <- merge(d13d18, sampledf, by=c('sample_id'), all=TRUE)

ent <- merge(sample_ent, entity, by=c('entity_id'), all=TRUE)

sitedf <- merge(ent, site, by=c('site_id'), all=TRUE)

df <- as.data.table(merge(sitedf, ages, by=c('sample_id'), all=TRUE))
df <- df[d18O == "yes"]


# ----------------------------------------------------------------------------- 
# SMOOTH AND GET STATISTICS
# ----------------------------------------------------------------------------- 
df.in <- df[ , c("d18O_measurement", "interp_age", "site_id", "latitude", "longitude")]
#
# get sites to loop across
these.ids <- unique(df.in$site_id)
# 
# function to rescale y-axis for plot checks
rescale_fun <- function(val, in_min, in_max, out_min, out_max){
  out_min + (val - in_min) * ((out_max - out_min) / (in_max - in_min))
}
#... initialize empty vectors 
myR <- mysd <- mysite <- meanAge <- vector()
idx <- 1
ageThreshold <- 0.8e3   #remove records with age ranges less than this (set to zero if you don't want this)
# --- run loop
for(i in 1:length(these.ids)){
  
  thisSite <- these.ids[i]
  #... grab the points from this location
  myPts <- df[which(df$site_id==thisSite), ]
  myPts <- myPts[order(myPts$interp_age, decreasing=TRUE), ]
  
  ageSpan <- max(na.omit(myPts$interp_age)) - min(na.omit(myPts$interp_age))
  
  if(all(is.na(myPts$interp_age)) | ageSpan < ageThreshold){ 
  } else{
    gapThreshold <- 2500
    myPts$gap <- c(NA, myPts[1:(length(myPts$interp_age)-1), ]$interp_age - myPts[2:length(myPts$interp_age), ]$interp_age)
    myPts$BeyondGap <- abs(myPts$gap) > gapThreshold
    gapDex <- which(myPts$BeyondGap==T)
    
    myPts <- as.data.table(cbind(myPts$interp_age, myPts$d18O_measurement, myPts$d13C_measurement,
                                 myPts$gap, myPts$BeyondGap))
    colnames(myPts) <- c('interp_age', 'd18O', 'd13C', 'gap', 'BeyondGap')
    
    ageSpan <- max(na.omit(myPts$interp_age)) - min(na.omit(myPts$interp_age))
    
    if(length(gapDex) < 1){
      
    } else{
      for(j in 1:length(gapDex)){
        #... collect the years 
        youngYear <- myPts[gapDex[j], ]$interp_age
        oldYear <- myPts[gapDex[j]-1, ]$interp_age
        theNAyears <- seq(oldYear-1, youngYear+1, by=-10)
        theAddition <- as.data.table(cbind(theNAyears, -9999, -9999, -9999, -9999))
        colnames(theAddition) <- c('interp_age', 'd18O', 'd13C', 'gap', 'BeyondGap')
        
        #... add to the data
        if(j==1){
          myPts_x <- rbind(myPts, theAddition)
        } else{
          myPts_x <- rbind(myPts_x, theAddition)
        }
      }
      
      myPts <- myPts_x
    }
    
    
    myAges <- seq(min(na.omit(myPts$interp_age)), max(na.omit(myPts$interp_age)), by=100)
    xRecord <- approxfun(x=myPts$interp_age, y=myPts$d18O, rule=1)
    xDF.temp <- as.data.table(cbind(myAges, xRecord(myAges))) ; colnames(xDF.temp) <- c('age', 'd18')
    xDF.temp$d18[which(xDF.temp$d18 < -100)] <- NA
    xDF.temp <- xDF.temp[order(xDF.temp$age, decreasing=TRUE), ]
    # get rooling mean (2000 year)
    xDF.list <- frollmean(xDF.temp, 20)
    xDF <- as.data.table(cbind(xDF.list[[1]], xDF.list[[2]])); colnames(xDF) <- c('age', 'd18')
    # --- see if it worked
    # insol_new <- rescale_fun(tmp.df[d18O>-100]$insol_wm2, min(na.omit(tmp.df[d18O>-100]$insol_wm2)),
    #                      max(na.omit(tmp.df[d18O>-100]$insol_wm2)), min(na.omit(tmp.df[d18O>-100]$d18O)),
    #                      max(na.omit(tmp.df[d18O>-100]$d18O)))
    # 
    # plot(xDF.temp$age, xDF.temp$d18, type='l', col='red')
    # lines(tmp.df[d18O>-100]$interp_age, ifelse((insol_new-mean(insol_new) >= 0), mean(insol_new) - abs(insol_new-mean(insol_new)),
    #                                            mean(insol_new) + abs(insol_new-mean(insol_new))), col='blue')
    # lines(xDF$age, xDF$d18, col='black')
    # title(paste("this site: ", thisSite, sep=''))
    
    # ... loop through 21 kyr intervals
    xDF.int <- xDF[!is.na(d18)]
    int_yrs <- 21e3
    n_int <- ceiling(ageSpan / int_yrs)
    int_max <- max(na.omit(xDF.int$age))
    int_min <- int_max - int_yrs
    keep_going <- TRUE
    while(keep_going==TRUE){
      tDF <- xDF.int[age < int_max & age >= int_min]
      # if data span > 8kyr, take stats
      if((max(tDF$age) - min(tDF$age)) > 8e3 & nrow(tDF > 40)){
        myR[idx] <- max(na.omit(tDF$d18)) - min(na.omit(tDF$d18))
        mysd[idx] <- sd(na.omit(tDF$d18))
        mysite[idx] <- thisSite
        meanAge[idx] <- mean(na.omit(tDF$age))
        # update idx
        idx <- idx + 1
      }
      # update int_max and min
      int_max <- max(xDF.int[age <= int_min]$age)
      int_min <- int_max - int_yrs
      if(int_max < min(na.omit(xDF.int$age)) + 8e3){
        keep_going <- FALSE
      }
    }
    
    
  }
  
  
}

dftest <- as.data.table(cbind(myR, mysd, mysite, meanAge)) ; colnames(dftest) <- c('xrange', 'xsd', 'site_id', 'meanAge')


dftest$percentile_range <- ecdf(dftest$xrange)(dftest$xrange)

# merge with lat long data
dfm <- merge(dftest, site, by='site_id')


# ------------------------------------------------------------------------
# --- SAVE RESULT 
saveDir <- getwd()
saveRDS(dfm, 
        file = paste(saveDir, 
                     "SISALv2_2kyrSmoothRangesINSOL_perOrbCycle.RDS", 
                     sep='/'))
# ------------------------------------------------------------------------









# --- READ IN DATA
setwd(saveDir)

# ... read in dat
df <- readRDS('SISALv2_2kyrSmoothRangesINSOL_perOrbCycle.RDS')

# WHERE TO SAVE
# save.here <- 'C:/Users/tkukl/OneDrive/Documents/Talks/AGU/2023/presentation'
save.here <- 'P:/iEBM_RUNS/SPATIAL_SCALES_PAPER/figure_panels/F1_spaceTimeScales'

# --- PLOT
# ... median line
med <- median(df[latitude > -37 & latitude < 37]$xrange)

# [1] HISTOGRAM
p.hist <- ggplot(df[latitude > -37 & latitude < 37]) +
  coord_cartesian(ylim=c(0,0.5)) +
  # median line
  annotate(geom='segment', x=med, xend=med, y=0, yend=0.5,
           size=1, color='darkred', linetype='dashed') +
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
ggsave(save.name, p.hist, width=18, height=15, units='cm')

# [2] EXAMPLE RECORD
brks = c(-15:-11)
coeff = -30
int_corr <- 7
min_age <- 194e3/1e3
max_age <- 245e3/1e3
p.example <- ggplot() +
  coord_cartesian(xlim=c(max_age, min_age), ylim=c(-9,-5.5)) + 
  geom_line(data=tmp.df[interp_age > min_age], aes(x=interp_age/1e3, y=((insol_wm2/coeff)+ int_corr)), 
            size=2, color='gray') +
  geom_line(data=xDF.temp[age > min_age], aes(x=age/1e3, y=d18), color='#9FC131', size=2) + 
  geom_line(data=xDF[age > min_age], aes(x=age/1e3, y=d18), color='#042940', size=3.5) +
  scale_y_continuous(
    # Features of the first axis
    name = expression("δ"^"18"*"O"~"(‰ VPDB)"),
    # Add a second axis and specify its features (adjust for coeff at labels, couldn't figure it out otherwise)
    sec.axis = sec_axis((~. - int_corr) , name=expression("DJF Insolation (W m"^"-2"*")"), breaks=brks,
                        labels=brks*coeff)
  ) +
  scale_x_continuous(name="Age (kyr bp)", expand=c(0,0)) + 
  theme_few() +
  theme(axis.text = element_text(size=20), axis.title=element_text(size=24), 
        axis.title.y.right = element_text(color = "gray"), axis.text.y.right = element_text(color="gray"))

save.fn <- "Cheng13_example.png"
save.name <- paste(save.here, save.fn, sep='/')
ggsave(save.name, p.example, width=20, height=13, units='cm')


