fig7 <- function(km_tsresults, km_pmrestults, csz) {
  cairo_ps(filename = paste(spath,'images/fig7.eps', sep=''))
  par(mar = c(5,5,2,5)) # Leave room on right side axis
  plot(km_tsresults$k,km_tsresults$wss2tss, type='b', lty=1, pch=19,
       col='black', cex=csz, cex.lab=csz, cex.axis=csz,
       axes=F, xlab=NA, ylab=NA, ylim=c(0,0.8))
  box()
  lines(km_pmresults$k,km_pmresults$wss2tss, type='b', lty=1, pch=2,
        col='black', cex=csz)
  axis(side = 1, at=K, label=K, col.axis='black', cex.axis=csz) # Draw x axis
  mtext(side = 1, line = 3, cex=csz, col='black', 'Number of Clusters')
  axis(side = 2, col.axis='black', cex.axis=csz) # Draw y1 axis
  mtext(side = 2, line = 3, cex=csz, col='black',
        'Cluster Dissimilarity [WSS/TSS]')
  par(new = T)
  # Plot compute time results
  y2max=max(c(km_tsresults$compute_secs,km_pmresults$compute_secs))
  plot(km_tsresults$k,km_tsresults$compute_secs, type='b', lty=2, pch=19,
       col='dodgerblue3',
       cex.axis=csz, cex=csz,
       axes=F, xlab=NA, ylab=NA, ylim=c(0, y2max))
  lines(km_pmresults$k,km_pmresults$compute_secs, type='b', lty=2, pch=2,
        col='dodgerblue3', cex=csz)
  y2lab=seq(0, y2max, length.out=10)
  y2lab=round(y2lab, digits=-1)
  axis(side = 4, at=y2lab, label=y2lab, col.axis='dodgerblue3', cex.axis=csz)
  mtext(side = 4, line = 3, cex=csz, col='dodgerblue3',
        'Compute Time [seconds]')
  legend("top", ncol=2, cex=csz*0.9, pt.cex=csz*0.9,
         legend=c('TS Sep.','PM Sep.',
                  'TS time', 'PM time'),
         lty=c(1,1,2,2), lwd=csz, pch=c(19,2,19,2),
         col=c('black', 'black', 'dodgerblue3', 'dodgerblue3'))
  dev.off()
  cat('Done\n')
}
