
cos.spheric.distance <- function(xy.pt, XY){
	# xy.pt: vector c(x, y)
	# XY: matrix X Y
	sin.y.p <- sin(xy.pt[2] * pi/180)
	sin.y.m <- sin(XY[, 2] * pi/180)
	cos.y.p <- cos(xy.pt[2] * pi/180)
	cos.y.m <- cos(XY[, 2] * pi/180)
	d.cos.x.pm <- cos((xy.pt[1] - XY[, 1]) * pi/180)
	sin.y <- sin.y.p * sin.y.m
	cos.y <- cos.y.p * cos.y.m
	cosd <- as.numeric(sin.y + cos.y * d.cos.x.pm)
	cosd[cosd > 1.0] <- 1.0
	return(cosd)
}

spheric.distance <- function(xy.pt, XY){
	# xy.pt: vector c(x, y)
	# XY: matrix X Y
	acos(cos.spheric.distance(xy.pt, XY))
}

cartesian.distance <- function(xy.pt, XY){
	# xy.pt: vector c(x, y)
	# XY: matrix X Y
	as.numeric(sqrt((xy.pt[1] - XY[, 1])^2 + (xy.pt[2] - XY[, 2])^2))
}

distance.vector <- function(xy.pt, XY, spheric){
	# xy.pt: vector c(x, y)
	# XY: matrix X Y
	if(spheric)
		spheric.distance(xy.pt, XY)
	else
		cartesian.distance(xy.pt, XY)
}

distance.matrix <- function(xy.pts, XY, spheric){
	# xy.pts: matrix X Y
	# XY: matrix X Y
	dst <- matrix(NA, nrow = nrow(XY), ncol = nrow(xy.pts))
	for (i in 1:nrow(xy.pts))
		dst[, i] <- distance.vector(as.numeric(xy.pts[i, ]), XY, spheric)
	return(dst)
}

.distance.pixel <- function(xy, XY, nmin, nmax, spheric){
	# xy: vector
	# XY: matrix (X, Y)
	# nmin, nmax: integer
	xdst <- distance.vector(xy, XY, spheric)
	odst <- order(xdst)
	xdst <- xdst[odst]
	nmean <- 1:as.integer((nmin + nmax)/2)
	rs.const <- mean(xdst[nmean])
	npts <- sum(xdst <= rs.const)
	rmax <- ifelse(npts < nmin, xdst[nmin], ifelse(npts > nmax, xdst[nmax], rs.const))
	npts <- sum(xdst <= rmax)

	list(dst = xdst, rmax = rmax, npts = npts, order = odst)
}
