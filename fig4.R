fig4 <- function(input, spy, lb, ub, csz) {
  ### input: Input vector or 1 column array of values (e.g. air temp)
  ### begyr: Calendar year of the first year of input data (e.g., 1995)
  ### spy: Number of samples (in this case weeks) per year
  library(PolarMetrics)
  dpy=365
  c=4 # Number of years of data
  input=input[1:208,]
  t <- as.numeric(strftime(index(input), format = '%j')) # Days of year
  r <- t2rad(t,dpy)     # Transform days of year to radians
  v <- as.vector(input) # Scalar values to be transformed
  vx <- mean(vec.x(r,v), na.rm=TRUE) # Avg horizontal vector
  vy <- mean(vec.y(r,v), na.rm=TRUE) # Avg vertical vector
  rvec <- vec_ang(vx,vy) # Angle of the avg vector (annual median)
  avec <- avec_ang(rvec)   # Vert opposite of med
  offsetDOY <- rad2d(avec,dpy) # Cal dy (1-365 beg at Jan 1) marking min
  offsetIdx <- rad2idx(avec,spy) # Index (1-spy) marking min
  AnnCum <- sum_cycle(v,offsetIdx,spy)$cumsum # Cum sum within each pheno yr
  # Note, output is re-centered, and accumulates starting from the annual
  # minimum (offsetIdx) and has c-1 yrs of data due to re-centering.
  c.AnnCum <- c-1       # Number of years in AnnCum
  cy <- 2 # Year within which to examine and plot cumulative values
  es <- window_idx(AnnCum,c.AnnCum,cy,lb,ub)[1] # Indx of AnnCum marking ES
  ms <- window_idx(AnnCum,c.AnnCum,cy,lb,ub)[3] # Indx of AnnCum marking MS
  ls <- window_idx(AnnCum,c.AnnCum,cy,lb,ub)[5] # Indx of AnnCum marking LS

  t2 <- seq(offsetDOY,by=7,length.out=spy+1)
  y <- c(0,AnnCum[(spy+1):(2*spy)])
  plot(t2,y,type='l', xlim=c(15,400),
    xaxt='n', ylim=c(0,max(y)),yaxs='i', yaxt='n', xlab='', ylab='',
    col='gray45', lwd=csz*2,
                cex=csz, cex.axis=csz, cex.lab=csz, cex.main=csz)
  y1 <- seq(0,30,by=5)
  y2 <- seq(15,85,by=10)
  x1 <- seq(0,360,by=30)
  x2 <- seq(0,360,by=30)
  axis(side=2, labels=y1, at=y1, padj=0.7, col='black', col.axis='black',
    cex.axis=csz) # Plot left y axis
  mtext('Cumulative Sum', col='black',
    side=2, line=1.9, cex=csz*0.9)
  ## x1 axis
  axis(side=1, labels=x1, at=x1+offsetDOY, col='black', col.axis='black',
    cex.axis=csz) # Plot right y axis
  ## x2 axis
  par(tcl = -0.5)       # Switch back to outward facing tick marks
  axis(side=3, labels=x2, at=x2, col='gray45', col.axis='gray45',
    cex.axis=csz) # Plot right y axis
  par(tcl = 1)          # Switch to inward facing tick marks
  mtext('Day of Calendar Year',side=3,line=-1.5,
    col='gray45',cex=csz*0.9)
  mtext('Day of Vector Centered Year',side=1,line=-1,
    col='black',cex=csz*0.9)
  ##
  axis(side=4, labels=y2, at=y2/(100/max(y)), col='gray45',
    col.axis='gray45',
    cex.axis=csz) # Plot right y axis
  mtext('% Cumulative Sum',side=4,line=-2,
    col='gray45',cex=csz*0.9)
  lines(t2[(es-spy+1):(ls-spy+1)],
    y[(es-spy+1):(ls-spy+1)], col='black', lwd=csz*3)
  #
  segments(offsetDOY,0,offsetDOY,0.95*max(y),
    col='gray45',lwd=csz*3, lty=2) # Phen yr offset dahsed line
  arrows(offsetDOY,0.95*max(y),offsetDOY,0.98*max(y),
    col='gray45',lwd=csz*3) # Phen yr offset arrow head
  text(offsetDOY+13,0.45*max(y),srt=90,
    'Anti-vector/Offset from Calendar Yr',cex=csz,col='brown3')
  #
  segments(t2[es-spy+1],y[(es-spy+1)],t2[es-spy+1],0.85*max(y),
    col='gray45',lwd=csz*3, lty=2) # ES vertical dashed line
  arrows(t2[es-spy+1],0.85*max(y),t2[es-spy+1],0.93*max(y),
    col='gray45',lwd=csz*3)        # ES arrow head (calendar yr)
  segments(t2[es-spy+1],y[(es-spy+1)],t2[es-spy+1],0.15*max(y),
    col='black',lwd=csz*3, lty=2)  # ES vertical dashed line
  arrows(t2[es-spy+1],0.15*max(y),t2[es-spy+1],0.07*max(y),
    col='black',lwd=csz*3)         # ES arrow head (phen yr)
  segments(t2[es-spy+1],y[(es-spy+1)],max(t2),y[(es-spy+1)],
    col='gray45',lwd=csz*3,lty=2)  # ES horizontal dashed line
  text(t2[es-spy+1]+30,y[(es-spy+1)]-1,srt=0,'ES',cex=csz,col='brown3')
  #
  segments(t2[ms-spy+1],y[(ms-spy+1)],t2[ms-spy+1],0.85*max(y),
    col='gray45',lwd=csz*3, lty=2) # MS vertical dashed line
  arrows(t2[ms-spy+1],0.85*max(y),t2[ms-spy+1],0.93*max(y),
    col='gray45',lwd=csz*3)        # MS arrow head (calendar yr)
  segments(t2[ms-spy+1],y[(ms-spy+1)],t2[ms-spy+1],0.15*max(y),
    col='black',lwd=csz*3, lty=2)  # MS vertical dashed line
  arrows(t2[ms-spy+1],0.15*max(y),t2[ms-spy+1],0.07*max(y),
    col='black',lwd=csz*3)         # MS arrow head (phen yr)
  segments(t2[ms-spy+1],y[(ms-spy+1)],max(t2),y[(ms-spy+1)],
    col='gray45',lwd=csz*3,lty=2)  # MS horizontal dashed line
  text(t2[ms-spy+1]+30,y[(ms-spy+1)]-1,srt=0,'MS',cex=csz,col='brown3')
  #
  segments(t2[ls-spy+1],y[(ls-spy+1)],t2[ls-spy+1],0.85*max(y),
    col='gray45',lwd=csz*3, lty=2) # LS vertical dashed line
  arrows(t2[ls-spy+1],0.85*max(y),t2[ls-spy+1],0.93*max(y),
    col='gray45',lwd=csz*3)        # LS arrow head (calendar yr)
  segments(t2[ls-spy+1],y[(ls-spy+1)],t2[ls-spy+1],0.15*max(y),
    col='black',lwd=csz*3, lty=2)  # LS vertical dashed line
  arrows(t2[ls-spy+1],0.15*max(y),t2[ls-spy+1],0.07*max(y),
    col='black',lwd=csz*3)         # LS arrow head (phen yr)
  segments(t2[ls-spy+1],y[(ls-spy+1)],max(t2),y[(ls-spy+1)],
    col='gray45',lwd=csz*3,lty=2)  # LS horiztontal dashed line
  text(t2[ls-spy+1]+30,y[(ls-spy+1)]-1,srt=0,'LS',cex=csz,col='brown3')
  #
  par(tcl = -0.5) # Switch back to outward facing tick marks
}
