sliding_window <- function(input, spy, breaks) {
  ### input: Input vector or 1 column array of values (e.g. air temp)
  ### spy: Number of samples (in this case weeks) per year
  library(PolarMetrics)
  vname=colnames(input)   # Variable name
  dpy <- 365              # Days/year
  t <- as.numeric(strftime(index(input), format = '%j')) # Days of year
  r <- t2rad(t,dpy)       # Transform days of year to radians
  v <- as.vector(input)
  lv <- length(v)
  VX <- vec.x(r,v)        # Horizontal vectors
  VY <- vec.y(r,v)        # Vertical vectors
  sqn=floor(seq(1, length(v), length.out=breaks))
  output=data.frame(fyr=rep(NA,breaks-1), Vavg=rep(NA,breaks-1),
                    Smag=rep(NA,breaks-1))
  for (I in 1:(length(sqn)-1)) {
    beg=sqn[I] # Start index of moving window
    end=beg+spy-1 # End index of 1-yr moving window
    vx <- mean(VX[beg:end], na.rm=T) # Avg horizontal vector
    vy <- mean(VY[beg:end], na.rm=T) # Avg vertical vector
    output$fyr[I] <- sqn[I]/spy
    output$Vavg[I] <- rad2d(vec_ang(vx,vy), dpy) # DOY of avg vec (50 pctile)
    output$Smag[I] <- vec_mag(vx,vy)     # Magnitude (length) of avg vector
  }
  return(output)
}
