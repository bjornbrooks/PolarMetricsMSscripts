fig1 <- function(input1, input2, csz) {
  par(xpd=TRUE) # Allow for drawing outside plot area
  input1=input1[1:208,]
  input2=input2[1:208,]
  t <- as.numeric(strftime(index(input1), format = '%j')) # Days of year
  t[53:104]=t[53:104]+365
  t[105:156]=t[105:156]+2*365
  t[157:208]=t[157:208]+3*365
  clrs <- c('blue3', 'red') # Point colors
  # Plot 1st var
  plot(t, as.vector(input1), type='l', lwd=csz,
       xlab='Days', ylab='', ylim=range(c(input1,input2)),
       cex=csz, cex.lab=csz, cex.axis=csz)
  mtext('Standardized Magnitude', side=2, line=2.5, cex=csz)
  lines(t, as.vector(input2), lwd=csz)
  points(t, as.vector(input1), pch=19, col=clrs[1],
         cex=csz)
  points(t, as.vector(input2), pch=19, col=clrs[2],
         cex=csz)
  legend('top', inset=c(0,-0.1), c('Temp','NEE'), ncol=2,
         col=clrs, lty=1, lwd=csz, pch=19, cex=csz)
}
