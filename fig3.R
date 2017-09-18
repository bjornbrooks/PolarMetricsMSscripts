fig3 <- function(input, spy, nrad_tick, csz) {
  library(plotrix)        # Load plotrix for plotting
  library(PolarMetrics)
  dpy=365
  input=input[1:208,]
  t <- as.numeric(strftime(index(input), format = '%j')) # Days of year
  r <- t2rad(t,dpy)       # Transform days of year to radians
  v=as.vector(input)
  VX <- vec.x(r,v)        # Horizontal vectors
  VY <- vec.y(r,v)        # Vertical vectors
  vx <- mean(VX, na.rm=T) # Avg horizontal vector
  vy <- mean(VY, na.rm=T) # Avg vertical vector
  vm <- vec_mag(vx,vy)     # Magnitude (length) of average vector
  rvec <- vec_ang(vx,vy)  # Angle of the avg vector (annual median)
  avec <- avec_ang(rvec) # Vert opposite of med (avg min/pheno yr start)
  # Plot
  lab.pos <- c(seq(from=0, to=2*pi-(2*pi/8), by=(2*pi)/8))[-3]
  olabs <- c('0', '1/4 pi', '1/2 pi', '3/4 pi', 'pi',
    '5/4 pi', '3/2 pi', '7/4 pi')[-3]
  s.pos <- pi/2 # Radial position to start plotting from
  par(cex.axis=csz, cex.lab=csz) # Set size of both axis labels
  # Plot 1st var
  radial.plot(as.vector(input), r,
    clockwise=TRUE, start=s.pos, label.pos=lab.pos, labels=olabs,
    rp.type='p', line.col='black', show.radial.grid=FALSE,
    grid.col='black', grid.unit='', cex=csz,
    radial.lim=pretty(c(0,1), n=nrad_tick))
  radial.plot(as.vector(input), r,
    clockwise=TRUE, start=s.pos,
    rp.type='s', point.symbols=20,
    point.col='blue3', cex=csz,
    radial.lim=pretty(c(0,1), n=nrad_tick),
    radlab=TRUE, add=TRUE) # Add colored dots
  legend("topleft", c('Temp'), col='blue3',
    lty=1, lwd=csz, pch=19, cex=csz, pt.cex=csz, box.lwd=0)
  radial.plot(c(0,1),c(0,avec),
    clockwise=TRUE,start=s.pos,rp.type='r',
    lwd=csz*1.5, line.col='mediumseagreen', add=TRUE) # avec, opposite angle of rvec
  radial.plot(c(0,1),c(0,rvec),
    clockwise=TRUE,start=s.pos,rp.type='r',
    lwd=csz*1.5,line.col='mediumseagreen',add=TRUE) # rvec, Angle of avg vec
  radial.plot(c(0,vm),c(rvec,rvec),
    clockwise=TRUE,start=s.pos,rp.type='r',
    lwd=csz*1.5,line.col='red',add=TRUE) # vm, vector magnitude
  radial.plot(c(0,vm),c(rvec,rvec),
    clockwise=TRUE,start=s.pos,rp.type='s',
    point.symbols='*',cex=csz*1.5,
    point.col='red',add=TRUE)   # vm, endpoints
  text(sin(avec)*1.15, cos(avec)*1.1,
    'AV', col='mediumseagreen', cex=csz) # Add text label
  text(sin(rvec)*1.2, cos(rvec)*1.1,
    'RV', col='mediumseagreen', cex=csz) # Add text label
  text(-0.2,0,
    'SMag',col='red',cex=csz) # Add text label
}
