get_modis <- function(dpath, cntr_ll, band, sp_aggr) {
  library(data.table)
  if (length(list.files(dpath, pattern='.asc$')) == 0) {
    library(MODISTools)
    library(xts)
    km_v=2; km_h=2 # Set span of domain (in km)
    span_hv <- c(km_v, km_h) # No. of km in x,y d imensions to acquire
    product <- c('MOD13Q1')
    band_name <- paste('250m_16_days_', band, sep='') # Data to query
    period <- data.frame(lat=cntr_ll[1], long=cntr_ll[2],
      start.date=2000, end.date=2016, id=1) # Time period
    # Download MODIS data
    MODISSubsets(LoadDat = period, Products = product, Bands = band_name,
                 Size = span_hv, SaveDir = dpath, StartDate = T)
  }
  downloaded.file <- list.files(path = dpath, pattern = '.asc$',
                                full.names = TRUE)
  modis_data=data.table::fread(downloaded.file, header = F) 

  nc=ncol(modis_data)
  if (sp_aggr == TRUE) { # Aggregate across all pixels
    output=xts(rowMeans(modis_data[,11:nc])/10^4, # Take spatial average
               as.POSIXct(gsub('A', '', modis_data$V8), format='%Y%j'),
               tz='America/Chicago')
  } else if (sp_aggr == FALSE) { # Keep pixels separate
    output=xts(modis_data[,11:nc]/10^4,
               as.POSIXct(gsub('A', '', modis_data$V8), format='%Y%j'),
               tz='America/Chicago')
  }

  return(output)
}
