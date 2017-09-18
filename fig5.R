fig5 <- function(input, spy, lcut, hcut, yr_type, begyr, yr, legend, nrad_tick, csz) {
  ### input: Input vector or 1 column array of values (e.g. air temp)
  ### begyr: Calendar year of the first year of input data (e.g., 1995)
  ### yr: The sequential year (1, 2, ...) within input to plot
  ### spy: Number of samples (in this case weeks) per year
  library(plotrix)        # Load plotrix for plotting
  library(PolarMetrics)
  clrVec=c('palegreen2','mediumseagreen','seagreen4')
  dpy <- 365              # Days/year
  vname=colnames(input)   # Variable name
  nyr <- nrow(input)/spy  # No. of years in input
  t <- as.numeric(strftime(index(input), format = '%j')) # Days of year
  r <- t2rad(t,dpy)       # Transform days of year to radians
  v <- as.vector(input)
  VX <- vec.x(r,v)        # Horizontal vectors
  VY <- vec.y(r,v)        # Vertical vectors
  vx <- mean(VX, na.rm=T) # Avg horizontal vector
  vy <- mean(VY, na.rm=T) # Avg vertical vector
  vm <- vec_mag(vx,vy)     # Magnitude (length) of average vector
  rvec <- vec_ang(vx,vy)  # Angle of the avg vector (polar median)
  avec <- avec_ang(rvec) # Vert opposite of med (avg min/polar yr start)
  offset_doy <- rad2d(avec,dpy) # Cal dy (1-365 beg at Jan 1) marking NDVI min
  offset_idx <- rad2idx(avec,spy) # Index (1-spy) marking NDVI min
  ann_cum <- sum_cycle(v,offset_idx,spy)$cumsum # Cum sums within each polar yr
  npy <- length(ann_cum)/spy # No. of complete years in output
  # Note, output is re-centered, and accumulates starting from
  # the annual NDVI minimum (PyStartIdx) and has nyr-1 yrs of
  # data due to re-centering.
  # early-, mid-, late-S ann_cum indices
  es_idx <- window_idx(ann_cum,npy,yr,lcut,hcut)[1] # Idx in ann_cum of ES milestone
  ms_idx <- window_idx(ann_cum,npy,yr,lcut,hcut)[3] # Idx in ann_cum of MS milestone
  ls_idx <- window_idx(ann_cum,npy,yr,lcut,hcut)[5] # Idx in ann_cum of LS milestone
  if (yr_type == 'cal_yr') {
    es_t <- t[es_idx + offset_idx - 1] # Julian day marking erly ssn
    ms_t <- t[ms_idx + offset_idx - 1] # Julian day marking mid ssn
    ls_t <- t[ls_idx + offset_idx - 1] # Julian day marking late ssn
  } else if (yr_type == 'offset_yr') {
    es_t <- t[es_idx] # Polar year day in t.s. marking erly ssn
    ms_t <- t[ms_idx] # Polar year day marking mid ssn
    ls_t <- t[ls_idx] # Polar year day marking late ssn
  }
  Smag <- vec_mag(mean(VX[es_idx:ls_idx], na.rm=TRUE),
    mean(VY[es_idx:ls_idx], na.rm=TRUE)) # Avg vec mag (length)
  # Plot
  bidx.yr=((yr-1)*spy)
  v.yr=v[bidx.yr:(bidx.yr+spy-1)] # magnitude polar coordinates for this year
  r.yr=r[bidx.yr:(bidx.yr+spy-1)] # angular polar coordinates for this year
  if (yr_type == 'cal_yr') {
    es.r=r[es_idx + offset_idx - 1] # Cal yr radial angle for early season milestone
    ms.r=r[ms_idx + offset_idx - 1] # Cal yr radial angle for early season milestone
    ls.r=r[ls_idx + offset_idx - 1] # Cal yr radial angle for early season milestone
  } else if (yr_type == 'offset_yr') {
    es.r=r[es_idx] # Polar yr radial angle for early season milestone
    ms.r=r[ms_idx] # Polar yr radial angle for early season milestone
    ls.r=r[ls_idx] # Polar yr radial angle for early season milestone
  }
  # Polar plot phenology variables
  s.pos <- pi/2 # Radial position to start plotting from
  lab.pos <- c(seq(from=0, to=2*pi-(2*pi/12), by=(2*pi)/12))[-4]
  olabs <- c(month.abb[seq(from=1, to=12)])[-4]
  par(cex.lab=csz*0.6) # Set size of labels of inner rings
  par(cex.axis=csz*0.4) # Set size of labels circling the outer
  radial.plot(v.yr, r.yr,
    clockwise=TRUE, start=s.pos, label.pos=lab.pos, labels=olabs,
    rp.type='p', line.col='black', lwd=csz,
    show.radial.grid=FALSE,
    grid.col='black', grid.unit='',
    radial.lim=pretty(range(v,na.rm=T), n=nrad_tick))
  radial.plot(v.yr, r.yr,
    clockwise=TRUE, start=s.pos,rp.type='s', point.symbols=20,
    point.col='black', cex=csz,
    radial.lim=pretty(range(v,na.rm=T), n=nrad_tick),
    radlab=TRUE, add=TRUE) # Add colored dots
  radial.plot(c(0,1), c(0,es.r),
    clockwise=TRUE, start=s.pos, rp.type='r',
    lwd=csz*2, line.col=clrVec[1], add=TRUE) # Avg ES angle
  radial.plot(c(0,1),c(0,ms.r),
    clockwise=TRUE,start=s.pos,rp.type='r',
    lwd=csz*2,line.col=clrVec[2],add=TRUE) # Angle of mean vec
  radial.plot(c(0,1),c(0,ls.r),
    clockwise=TRUE,start=s.pos,rp.type='r',
    lwd=csz*2,line.col=clrVec[3],add=TRUE) # Avg LS angle
  radial.plot(c(0,1),c(0,avec),
    clockwise=TRUE,start=s.pos,rp.type='s',
    point.symbols='*',cex=csz,
    add=TRUE) # avec, opposite angle of rvec
  radial.plot(c(0,Smag),c(ms.r,ms.r),
    clockwise=TRUE,start=s.pos,rp.type='r',
    lwd=csz,line.col='red',add=TRUE) # vm, vector magnitude
  radial.plot(c(0,Smag),c(ms.r,ms.r),
    clockwise=TRUE,start=s.pos,rp.type='s',
    point.symbols='*',cex=csz,
    point.col='red',add=TRUE)   # vm, Magnitude of avg vec
  text(sin(avec)*1.15, cos(avec)*1.1,
    'AV', col='black', cex=csz) # Add text label
  text(sin(rvec)*1.2, cos(rvec)*1.1,
       'RV', col='black', cex=csz) # Add text label
  text(sin(es.r)*1.12, cos(es.r)*1.1,
    'ES', col=clrVec[1], cex=csz) # Add text label
  text(sin(ms.r)*1.12, cos(ms.r)*1.1,
       'MS', col=clrVec[3], cex=csz) # Add text label
  text(sin(ls.r)*1.12, cos(ls.r)*1.1,
    'LS', col=clrVec[3], cex=csz) # Add text label
  text(0,0.1,'SMag',col='red',cex=csz) # Add text label
  if (vname=='NDVI') {snt='NDVI'
  } else if (vname=='NEE_PI') {snt='NEE'
  } else if (vname=='P') {snt='Prec'
  } else if (vname=='PPFD_IN') {snt='PPFD'
  } else if (vname=='TA') {snt='AirT'
  } else {snt=vname
  }
  text(sin(5.5)*1.3, cos(5.7)*1.3,
       snt,col='black',cex=csz*1.25) # Add var name
  if (legend == TRUE) {
    legend('bottomright', inset= -0.05,
      as.character(begyr+yr-1),
      col=c('red'), pch=20,cex=0.65,pt.cex=2,box.lwd=0)
  }
}
