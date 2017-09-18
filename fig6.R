fig6 <- function(input1, input2, csz) {
  # Plot
  cairo_ps(filename = paste(spath,'images/fig6.eps', sep=''))
  par(mar = c(5,5,2,5)) # Leave room on right side axis
  plot(input1$fyr,input1$Vavg, type='b', lty=1, pch=19,
    col='mediumseagreen', cex=csz, cex.lab=csz, cex.axis=csz,
    axes=F, xlab=NA, ylab=NA, ylim=c(80,310))
  box()
  lines(input2$fyr,input2$Vavg, type='b', lty=1, pch=2,
        col='mediumseagreen', cex=csz)
  xlab=seq(1995, 2015, by=2)
  axis(side = 1, at=xlab, label=xlab, col.axis='black',
       cex.axis=csz) # Draw x axis
  mtext(side = 1, line = 3, cex=csz, col='black', 'Year')
  axis(side = 2, col.axis='mediumseagreen', cex.axis=csz) # Draw y1 axis
  mtext(side = 2, line = 3, cex=csz, col='mediumseagreen', 'Timing [DOY]')
  par(new = T)
  plot(input1$fyr,input1$Smag, type='b', lty=2, pch=19, col='firebrick2',
    cex.axis=csz, cex=csz,
    axes=F, xlab=NA, ylab=NA, ylim=c(0,0.3))
  lines(input2$fyr,input2$Smag, type='b', lty=2, pch=2,
        col='firebrick2', cex=csz)
  y2lab=seq(0,0.2, by=0.1)
  axis(side = 4, at=y2lab, label=y2lab, col.axis='firebrick2', cex.axis=csz)
  mtext(side = 4, line = 3, adj=0.25, cex=csz, col='firebrick2',
        'Seasonality [0-1]')
  legend("top", ncol=2, cex=csz*0.9, pt.cex=csz*0.9,
         legend=c('NEE Timing', 'P Timing',
                  'NEE Ssn', 'P Ssn'),
         lty=c(1,1,2,2), lwd=csz, pch=c(19,2,19,2),
         col=c('mediumseagreen', 'mediumseagreen', 'firebrick2', 'firebrick2'))
  dev.off()
}
