fig2 <- function(input1, input2, spy, nrad_tick, csz) {
	### input: Input vector or 1 column array of values (e.g. air temp)
	### begyr: Calendar year of the first year of input data (e.g., 1995)
	### yr: The sequential year (1, 2, ...) within input to plot
	### spy: Number of samples (in this case weeks) per year
	library(plotrix)        # Load plotrix for plotting
	library(PolarMetrics)
	dpy=365
	input1=input1[1:208,]
	input2=input2[1:208,]
	t <- as.numeric(strftime(index(input1), format = '%j')) # Extract days of year
	r <- t2rad(t,dpy)       # Transform days of year to radians
	# Plot
	lab.pos <- c(seq(from=0, to=2*pi-(2*pi/8), by=(2*pi)/8))[-3]
	olabs <- c('0', '1/4 pi', '1/2 pi', '3/4 pi', 'pi',
		'5/4 pi', '3/2 pi', '7/4 pi')[-3]
	clrs <- c('blue3', 'red') # Point colors
	s.pos <- pi/2 # Radial position to start plotting from
	par(cex.axis=csz, cex.lab=csz) # Set size of both axis labels
	# Plot 1st var
	radial.plot(as.vector(input1), r,
		clockwise=TRUE, start=s.pos, label.pos=lab.pos, labels=olabs,
		rp.type='p', line.col='black', show.radial.grid=FALSE,
		grid.col='black', grid.unit='', cex=csz,
		radial.lim=pretty(c(0,1), n=nrad_tick))
	radial.plot(as.vector(input1), r,
		clockwise=TRUE, start=s.pos,
		rp.type='s', point.symbols=20,
		point.col=clrs[1], cex=csz,
		radial.lim=pretty(c(0,1), n=nrad_tick),
		radlab=TRUE, add=TRUE) # Add colored dots
	# Plot 2nd var
	radial.plot(as.vector(input2), r,
		clockwise=TRUE, start=s.pos,
		rp.type='p', line.col='black', show.radial.grid=FALSE,
		grid.col='black', grid.unit='',
		radial.lim=pretty(c(0,1), n=nrad_tick), add=TRUE)
	radial.plot(as.vector(input2), r,
		clockwise=TRUE, start=s.pos,rp.type='s', point.symbols=20,
		point.col=clrs[2], cex=csz,
		radial.lim=pretty(c(0,1), n=nrad_tick),
		radlab=TRUE, add=TRUE) # Add colored dots
	legend("topleft", c('Temp','NEE'), col=clrs,
	       lty=1, lwd=csz, pch=19, cex=csz, pt.cex=csz, box.lwd=0)
}
