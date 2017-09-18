fig7kmeans <- function(lat, lon, K) {
  ndvi2 <- get_modis(dpath2, c(lat, lon), band='NDVI', sp_aggr=FALSE)
  ndvi2[paste(begyr, '/', begyrMODIS-1, sep='')] <- 0 # Set bad data year to
                                                # zero to assist interpolation
  # Create time index
  tim <- seq.Date(from=as.Date(paste(2000,1,7,sep='-')),
                 to=as.Date(paste(2000,12,31,sep='-')),
                 by=7) # Make time index
  for (yr in 2001:2016) {
    ti <- seq.Date(from=as.Date(paste(yr,1,7,sep='-')),
                   to=as.Date(paste(yr,12,31,sep='-')),
                   by=7) # Make time index
    tim <- c(tim, ti)
  }
  tim <- as.POSIXct(tim, tz='America/Chicago') # Set time attributes
  ndvi2 <- na.spline(merge(ndvi2,tim)) # Down-scale NDVI through
                                  # interpolation & merge with lef.7dy
  ndvi2[paste(begyr, '/', begyrMODIS-1, sep='')] <- NA # Set bad data yr to NA
  ndvi2 <- ndvi2[tim[(wpy+1):length(tim)],] # Extract 2001-2016
  ndvi2[ndvi2 < 0] <- 0
  # Make time series testing data for clustering
  test_ts <- ndvi2[1:(nrow(ndvi2)-wpy),] # Remove 2016 to match metrics output
  test_ts <- matrix(as.matrix(test_ts),
                     nrow=wpy, byrow=F) # Reform into annual chunks
  test_ts <- aperm(test_ts, c(2,1)) # Rotate matrix so each week is a variable
  test_ts <- normaliz(test_ts) # Standardize globally to a range of 0-1
  # Create testing data for polar metrics technique
  for (I in 1:ncol(ndvi2)) {
    if (I == 1) {
      test_pm <- calc_metrics(ndvi2[, I], yr_type=tme_sc,
                               spc=wpy, lcut=lb, hcut=ub, return.vecs=FALSE)
    } else {
      tmp <- calc_metrics(ndvi2[, I], yr_type=tme_sc,
                               spc=wpy, lcut=lb, hcut=ub, return.vecs=FALSE)
      test_pm <- rbind(test_pm, tmp)
    }
  }
  test_pm <- test_pm[,-1] # Remove year variable
  for (I in 1:ncol(test_pm)) {
    test_pm[, I]<- normaliz(test_pm[, I]) # Standardize each metric to a range of 0-1
  }
  rm(tmp)
  km_pmresults <- data.frame(k = K, compute_secs = rep(NA,11),
                             wss2tss= rep(NA,11))
  for (I in 1:11) {
    k <- K[I]
    t0 <- proc.time()
    km  <- kmeans(test_pm, centers=k, nstart=10^3, iter.max=10^2)
    km_pmresults$compute_secs[I] <- (proc.time()-t0)[3]
    # Ratio of total variance within clusters to total variance in data set,
    #   i.e., average representativeness of each cluster
    km_pmresults$wss2tss[I] <- km$tot.withinss/km$totss
    print(km_pmresults[I,])
  }
  rm(t0,k,km)
  km_tsresults <- data.frame(k = K, compute_secs = rep(NA,11),
                             wss2tss= rep(NA,11))
  for (I in 1:11) {
    k <- K[I]
    t0 <- proc.time()
    km  <- kmeans(test_ts, centers=k, nstart=10^3, iter.max=10^2)
    km_tsresults$compute_secs[I] <- (proc.time()-t0)[3]
    # Ratio of total variance within clusters to total variance in data set,
    #   i.e., average representativeness of each cluster
    km_tsresults$wss2tss[I] <- km$tot.withinss/km$totss
    print(km_tsresults[I,])
  }
  rm(t0,k,km)
  save(km_pmresults, km_tsresults, file=file.path(spath,'fig7kmeans.RData'))
}
