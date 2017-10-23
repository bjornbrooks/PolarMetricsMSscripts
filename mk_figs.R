rm(list=ls())

fig7newdata <- FALSE
spath <- '.'
WLEFlon <- -90.2723 # Longitude of WLEF
WLEFlat <- 45.9459 # Latitude of WLEF
AMFname <- 'AMF_US-PFa_BASE_HR_9-1.csv' # Name of Ameriflux data file
refcols <- c('NEE_PI', 'TA', 'PPFD_IN', 'P') # Variables to extract
begyr <- 1995 # First year of data
endyr <- 2015 # Last year of data
begyrMODIS <- 2001 # First complete year of MODIS NDVI data
dpy <- 365 # Days per year
wpy <- 52 # Weeks per year
tme_sc <- 'cal_yr' # Timing metrics output as 'cal_yr' or 'offset_yr'
dpath1 <- paste(spath,'ameriflux/',sep='') # location of AmeriFlux data
dpath2 <- paste(spath,'modis/',sep='') # location to download MODIS data to
lb <- 0.15 # [0:0.5] Lower bound of early season (ES) timing metric
ub <- 0.85 # [0.5:1] Upper bound of late season (LS) timing metric

# Load libraries
library(xts)

# Load custom functions
source(file.path(spath,'get_modis.R'))

# Custom function for normalizing values between 0 and 1
normaliz <- function(x){
  (x - min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))
}

### Load and parse data
cat('Loading and parsing input...')
lef <- read.csv(file=file.path(dpath1, AMFname)) # Load tower data
lef <- lef[, c('TIMESTAMP_START',
               refcols)] # Make subset of time, NEE, TA, PPFD_IN, P
lef[lef == -9999] = NA
time.start <-  as.POSIXct(gsub('(.{8})', '\\1 ', lef$TIMESTAMP_START),
                         format='%Y%m%d %H%M',
                         tz='America/Chicago') # Make data an XTS object
miss <- which(is.na(time.start)) # Tag rows with corrupt time data
if (length(miss) > 0) { # If missing values found then remove first
  lef <- xts(lef[-miss, 2:ncol(lef)],
          time.start[-miss]) # Make XTS object from good time tags
} else {
  lef <- xts(lef[, 2:ncol(lef)], time.start) # Make XTS object
}
lef$NEE_PI <- -lef$NEE_PI # Reverse the sign of NEE so + means uptake of C
lef <- lef[! is.na(index(lef))] # Remove dates with bad time data if any
lef.dy <- apply.daily(lef, FUN=mean, na.rm=T) # Calculate daily averages

### Loop over all yrs and average over 7 day periods beginning at Jan 1
### until 52 7-day periods are reached in each year. Note that the remaining
### 2-4 days are discarded.
for (yr in begyr:endyr) {
  tmp <- lef.dy[paste(yr)]
  ep <- seq(from=1, to=365, by=7) # Define 7-day sampling points
  tmp2 <- period.apply(tmp, INDEX=ep, FUN=mean,
                       na.rm=T)[1:wpy,] # 7-day avgs
  ti <- seq.Date(from=as.Date(paste(yr,1,7,sep='-')),
                 to=as.Date(paste(yr,12,31,sep='-')),
                 by=7) # Make time index
  ti <- as.POSIXct(ti, tz='America/Chicago') # Set time attributes
  if (yr == begyr) {
    lef.7dy <- xts(tmp2, ti, tz='America/Chicago')
  } else {
    lef.7dy <- rbind(lef.7dy, xts(tmp2, ti))
  }
}
rm(yr, ep, ti, tmp, tmp2)

# Gap-fill weeks with missing values with the average for that week
for (I in 1:ncol(lef.7dy)) {
  tmp  <-  matrix(lef.7dy[,I], ncol=(endyr-begyr+1))
  miss.idx <- which(is.na(tmp)) # Locate missing vals
  fill <- rowMeans(tmp, na.rm=T) # Expected average for each week
  fill.idx <- miss.idx %% wpy # Indices corresponding to fill vals (means)
  fill.idx[fill.idx == 0] <- wpy
  tmp[miss.idx] <- fill[fill.idx] # Insert avg's where missvals were
  lef.7dy[,I] <- tmp
  rm(tmp) # Clean up
}
rm(fill, fill.idx) # Clean up

# Download and load MODIS data
ndvi <- get_modis(dpath2, c(WLEFlat, WLEFlon), band='NDVI', sp_aggr=TRUE)
ndvi[paste(begyr, '/', begyrMODIS-1, sep='')] <- 0 # Set bad data year to
                                              # zero to assist interpolation
ndvi <- na.spline(merge(ndvi,index(lef.7dy))) # Down-scale NDVI through
                                # interpolation & merge with lef.7dy
ndvi[paste(begyr, '/', begyrMODIS-1, sep='')] <- NA # Set bad data yr to NA
ndvi <- ndvi[index(lef.7dy),] # Extract only the dates corresponding to lef.7dy
lef.7dy$NDVI <- ndvi
rm(ndvi) # Clean up
lef.7dy <- lef.7dy[,c('NDVI', refcols)] # Re-order columns for plotting

# Standardize all environmental variables to a range of 0-1
#   (Just so we can compare across variables).
for (I in 1:ncol(lef.7dy)) {
  lef.7dy[,I] <- normaliz(lef.7dy[,I])
}
cat('Done\n')

### Figure 1
cat('Processing figure 1...')
# Plot Fig. 1
source(file.path(spath,'fig1.R'))
#png(filename = paste(spath,'images/fig1.png', sep=''),
#  width = 600, height = 400, units = 'px')
cairo_ps(filename = paste(spath,'images/fig1.eps', sep=''))
fig1(lef.7dy$TA, lef.7dy$NEE_PI, csz=1.8)
dev.off()
cat('Done\n')

### Figure 2
cat('Processing figure 2...')
# Plot Fig. 2
source(file.path(spath,'fig2.R'))
#png(filename = paste(spath,'images/fig2.png', sep=''),
#  width = 600, height = 600, units = 'px')
cairo_ps(filename = paste(spath,'images/fig2.eps', sep=''))
fig2(lef.7dy$TA, lef.7dy$NEE_PI, spy=wpy, nrad_tick=4, csz=1.75)
dev.off()
cat('Done\n')

### Figure 3
cat('Processing figure 3...')
# Plot Fig. 3
source(file.path(spath,'fig3.R'))
#png(filename = paste(spath,'images/fig3.png', sep=''),
#  width = 600, height = 600, units = 'px')
cairo_ps(filename = paste(spath,'images/fig3.eps', sep=''))
fig3(lef.7dy$TA, spy=wpy, nrad_tick=4, csz=1.75)
dev.off()
cat('Done\n')

### Figure 4
cat('Processing figure 4...')
# Plot Fig. 4
source(file.path(spath,'fig4.R'))
#png(filename = paste(spath,'images/fig4.png', sep=''),
#  width = 600, height = 600, units = 'px')
cairo_ps(filename = paste(spath,'images/fig4.eps', sep=''))
fig4(lef.7dy$TA, spy=wpy, lb=0.15, ub=0.85, csz=1.5)
dev.off()
cat('Done\n')

### Figure 5
cat('Processing figure 5...')
### Calculate annual polar metrics for all environmental variables
for (I in 1:ncol(lef.7dy)) {
  if (I ==1) {
    MODISyrs <- paste(begyrMODIS, '/', endyr, sep='')
    assign(colnames(lef.7dy)[I],
           calc_metrics(lef.7dy[MODISyrs, I], yr_type=tme_sc,
                       spc=wpy, lcut=lb, hcut=ub, return.vecs=FALSE))
    rm(MODISyrs) # Clean up
  } else {
    assign(colnames(lef.7dy)[I],
           calc_metrics(lef.7dy[, I], yr_type=tme_sc,
                       spc=wpy, lcut=lb, hcut=ub, return.vecs=FALSE))
  }
}

# Make Fig. 5 plots
dfr=data.frame(yr=TA$yr, NDVIes=c(rep(NA,6),NDVI$es),
               NEE_PIes=NEE_PI$es[1:20], TAes=TA$es,
               PPFDes=PPFD_IN$es, Pes=P$es) # Data for linear plot
clrs=c('green3', 'seagreen2', 'firebrick1', 'orange1', 'blue3')
# Make polar plots for all environmental variables 1 year at a time
csz=7
csz2=3.5
source(file.path(spath,'fig5.R'))
for (yr in 1:(length(begyr:endyr)-1)) {
  yr0=sprintf('%02i',yr)
  png(filename = paste(spath,'images/fig5_', yr0, '.png', sep=''),
      width = 1400, height = 2000, units = 'px')
  par(mfrow=c(3,2))
  for (I in 1:ncol(lef.7dy)) {
    if (yr < 7 & I == 1) {
      plot.new()
      text(0.45,0.6, 'No', cex=csz*1.1)
      text(0.45,0.45, 'NDVI', cex=csz*1.1)
    } else {
      fig5(lef.7dy[, I], spy=wpy, lcut=lb, hcut=ub, yr_type=tme_sc,
           begyr, yr, legend=FALSE, nrad_tick=3, csz=csz)
    }
  }
  # Add linear time series plot
  par(mar = c(5,5,2,5)) # Leave room on right side axis
  plot(dfr$yr[1:yr]+begyr-1, dfr$NDVIes[1:yr], type='b', lty=1, pch=1,
       col=clrs[1], cex=csz2, lwd=csz2, cex.axis=csz2,
       xlab=NA, ylab=NA, xlim=c(begyr,endyr), ylim=c(50,220))
  lines(dfr$yr[1:yr]+begyr-1, dfr$NEE_PIes[1:yr], col=clrs[2], type='b',
        cex=csz2, lwd=csz2, lty=2, pch=2)
  lines(dfr$yr[1:yr]+begyr-1, dfr$TAes[1:yr], col=clrs[3], type='b',
        cex=csz2, lwd=csz2, lty=3, pch=3)
  lines(dfr$yr[1:yr]+begyr-1, dfr$PPFDes[1:yr], col=clrs[4], type='b',
        cex=csz2, lwd=csz2, lty=4, pch=4)
  lines(dfr$yr[1:yr]+begyr-1, dfr$Pes[1:yr], col=clrs[5], type='b',
        cex=csz2, lwd=csz2, lty=5, pch=5)
  mtext(side = 1, line = -1.3, cex=csz2*0.8, 'Year')
  mtext(side = 4, line = 2.5, cex=csz2*0.8, 'ES [Day of Year]')
  legend("topleft", cex=csz2*0.55, pt.cex=csz2, lwd=csz2, ncol=2,
         legend=c('NDVI','NEE','AirT','PPFD','Pr'),
         lty=1:5, pch=1:5, col=clrs)
  dev.off()

}
rm(dfr, clrs, csz, yr, yr0)
cat('Done\n')


### Figure 6
cat('Processing figure 6...')
# Use a custom sliding window function to look for non stationarity
#   in the polar metrics.
source(file.path(spath,'sliding_window.R'))
P.sw=sliding_window(lef.7dy$P, spy=wpy, breaks=50)
P.sw$fyr=P.sw$fyr+1995 # Set to calendar years
P.sw=P.sw[P.sw$fyr<2015,]
NEE.sw=sliding_window(lef.7dy$NEE_PI, spy=wpy, breaks=50)
NEE.sw$fyr=NEE.sw$fyr+1995 # Set to calendar years
NEE.sw=NEE.sw[NEE.sw$fyr<2015,]
TA.sw=sliding_window(lef.7dy$TA, spy=wpy, breaks=50)
TA.sw$fyr=TA.sw$fyr+1995 # Set to calendar years
TA.sw=TA.sw[TA.sw$fyr<2015,]
NDVI.sw=sliding_window(lef.7dy$NDVI, spy=wpy, breaks=50)
NDVI.sw$fyr=NDVI.sw$fyr+1995 # Set to calendar years
NDVI.sw=NDVI.sw[NDVI.sw$fyr<2015,]
NDVI.sw[NDVI.sw$fyr<2001,2:3]=NA # Fix bad data

# Explore regional precipitation data from NARR
#source(file.path(spath,'get_narr.R'))
#narr <- get_narr(varname='prate.sfc',
#                 lon=WLEFlon, lat=WLEFlat, begyr=begyr, endyr=endyr)
load(file.path(spath,'narr/narr.RData'))
time <- as.POSIXct(narr$datetime,
                   format = "%Y_%m_%d_%H", tz='America/Chicago')
narr.xt <- xts(narr[,2], time)
colnames(narr.xt)[1] <- 'P'
narr.dy=apply.daily(narr.xt, FUN=mean, na.rm=T) # Calculate daily averages
### Loop over all yrs and average over 7 day periods beginning at Jan 1
### until 52 7-day periods are reached in each year. Note that the remaining
### 2-4 days are discarded.
for (yr in begyr:endyr) {
  tmp=narr.dy[paste(yr)]
  ep <- seq(from=1, to=365, by=7) # Define 7-day sampling points
  tmp2 <- period.apply(tmp, INDEX=ep, FUN=mean,
                       na.rm=T)[1:wpy,] # 7-day avgs
  ti <- seq.Date(from=as.Date(paste(yr,1,7,sep='-')),
                 to=as.Date(paste(yr,12,31,sep='-')),
                 by=7) # Make time index
  ti <- as.POSIXct(ti, tz='America/Chicago') # Set time attributes
  if (yr == begyr) {
    narr.7dy <- xts(tmp2, ti, tz='America/Chicago')
  } else {
    narr.7dy <- rbind(narr.7dy, xts(tmp2, ti))
  }
}
rm(yr, ep, tmp, tmp2)
colnames(narr.7dy)='P'
narr.7dy$P <- normaliz(narr.7dy$P) # Standardize to a range of 0-1
# Calculate annual polar metrics for all environmental variables
P.narr <- calc_metrics(narr.7dy$P, yr_type=tme_sc,
                      spc=wpy, lcut=lb, hcut=ub, return.vecs=FALSE)

# Sliding window
P.narr.sw=sliding_window(narr.7dy$P, spy=wpy, breaks=50)
P.narr.sw$fyr=P.narr.sw$fyr+1995 # Set to calendar years
P.narr.sw=P.narr.sw[P.narr.sw$fyr<2015,]

# Plot
source(file.path(spath,'fig6.R'))
fig6(input1= NEE.sw, input2=P.sw, csz=1.7)
cat('Done\n')

### Figure 7
K = 2^(1:11) # Number of cluster centers to try
if (fig7newdata == TRUE) {
  cat('Downloading and clustering new data for figure 7...')
  # Download MODIS and cluster modis data
  source(file.path(spath,'fig7kmeans.R'))
  fig7kmeans(lat=WLEFlat, lon=WLEFlon, K = K)
  cat('Done\n')
}
cat('Loading figure 7 data...')
load(file.path(spath,'fig7kmeans.RData'))
cat('Done\n')

cat('Processing figure 7...')
# Plot
source(file.path(spath,'fig7.R'))
fig7(km_tsresults, km_pmrestults, csz=1.7)
cat('Done\n')

print('All Done!')
