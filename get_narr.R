get_narr <- function(varname, lon, lat, begyr, endyr) {
  library(RNCEP)
  ## Retrieve the temperature from a particular pressure level for
  ## a specified spatial and temporal extent
  narr <- NCEP.gather(variable=varname, level='gaussian',
                        months.minmax=c(1,12), years.minmax=c(begyr,endyr),
                        lat.southnorth=c(lat,lat), lon.westeast=c(lon,lon),
                        reanalysis2 = FALSE, return.units = TRUE)
  ## Then convert the 3-D array to a data.frame ##
  df <- NCEP.array2df(wx.data=narr, var.names=varname)
  df.lons=unique(df$longitude)
  df.lats=unique(df$latitude)
  lon.nn=df.lons[which.min(abs(lon-df.lons))]
  lat.nn=df.lats[which.min(abs(lat-df.lats))]
  df <- df[which(df$latitude == lat.nn & df$longitude == lon.nn),]
  print(paste('Download coords: ', lon.nn, '(lon),', lat.nn,'(lat)'))
  output <- data.frame(datetime=df[,1], V=df[,4])

  return(output)
}