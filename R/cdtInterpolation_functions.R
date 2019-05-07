
#' CDT interpolation methods
#' 
#' This function interpolates stations data on a rectangular grid in the x-y plane.
#' The radius of influence for the interpolation method varies for each grid point. 
#' 
#' @param locations A matrix containing the coordinates of the station data. The 1st column is the longitude and the 2nd coulmn the latitude. 
#' @param values A vector containing the values of the station. The length of \code{values} must be equal to the number of row of \code{locations}.
#' @param newdata A matrix containing the coordinates of the grid to interpolate. The 1st column is the longitude and the 2nd coulmn the latitude.
#' @param nmin The minimum number of stations to be used to interpolate a grid point.
#' @param nmax The maximum number of stations to be used to interpolate a grid point.
#' @param method Method to calculate weights for the interpolation. Should be \strong{"idw"} (default), \strong{"shepard"}, \strong{"spheremap"} or \strong{"kriging"}.
#' \itemize{
#' 	\item \strong{"idw"}: Inverse distance weighted
#' 	\item \strong{"shepard"}: Modified Shepard interpolation
#' 	\item \strong{"spheremap"}: Spheremap interpolation method
#' 	\item \strong{"kriging"}: Ordiranry kriging
#' }
#' @param spheric If \code{FALSE} (default), then a cartesian distance will be computed. If set to \code{TRUE}, a spherical distance will be computed.
#' @param vgm A \code{gstat} variogram model. Example: \code{vgm(10.2, "Sph", 2.7, 1.1)}.
#' @param p The power to use in weight calculation.
#' 
#' @return A matrix containing the interpolated data. The 1st and 2nd columns are the same as \code{newdata} and the 3rd column contains the interpolated values.
#' 
#' @export

cdtInterp <- function(locations, values, newdata, nmin = 4, nmax = 10,
					method = "idw", spheric = FALSE, vgm = NULL, p = 2.0)
{
	if(nrow(locations) != length(values)) stop("The length of 'values' must be equal to the number of row of 'locations'")
	if(method == "kriging" && is.null(vgm)) stop("No variogram model found.")
	locations <- as.matrix(locations[, 1:2, drop = FALSE])
	newdata <- as.matrix(newdata[, 1:2, drop = FALSE])

	interp.res <- switch(method,
		"idw" = idw.interp(locations, values, newdata, nmin, nmax, spheric, p),
		"shepard" = shepard.interp(locations, values, newdata, nmin, nmax, spheric, p),
		"spheremap" = spheremap.interp(locations, values, newdata, nmin, nmax, spheric),
		"kriging" = kriging.interp(locations, values, newdata, vgm, nmin, nmax, spheric),
		cat("Unknown interpolation method\n")
	)

	return(interp.res)
}

##########################################################################################

idw.interp <- function(locations, values, newdata, nmin, nmax, spheric, p = 2){
	# locations: matrix X Y
	# values: vector
	# newdata: matrix X Y
	# nim: integer
	# nmax: integer
	# spheric: logical
	# p: real
	ngrd <- nrow(newdata)
	interp.res <- matrix(NA, nrow = ngrd, ncol = 3)
	interp.res[, 1:2] <- newdata
	for(j in seq(ngrd)){
		xy.g <- as.numeric(newdata[j, ])
		interp.res[j, 3] <- .idw.pixel(xy.g, locations, values, nmin, nmax, spheric, p)
	}

	return(interp.res)
}

shepard.interp <- function(locations, values, newdata, nmin, nmax, spheric, p = 2){
	ngrd <- nrow(newdata)
	interp.res <- matrix(NA, nrow = ngrd, ncol = 3)
	interp.res[, 1:2] <- newdata
	for(j in seq(ngrd)){
		xy.g <- as.numeric(newdata[j, ])
		interp.res[j, 3] <- .shepard.pixel(xy.g, locations, values, nmin, nmax, spheric, p)
	}

	return(interp.res)
}

barnes.interp <- function(locations, values, newdata, nmin, nmax, spheric, p = 2){
	ngrd <- nrow(newdata)
	interp.res <- matrix(NA, nrow = ngrd, ncol = 3)
	interp.res[, 1:2] <- newdata
	for(j in seq(ngrd)){
		xy.g <- as.numeric(newdata[j, ])
		interp.res[j, 3] <- .barnes.pixel(xy.g, locations, values, nmin, nmax, spheric, p)
	}

	return(interp.res)
}

cressman.interp <- function(locations, values, newdata, nmin, nmax, spheric){
	ngrd <- nrow(newdata)
	interp.res <- matrix(NA, nrow = ngrd, ncol = 3)
	interp.res[, 1:2] <- newdata
	for(j in seq(ngrd)){
		xy.g <- as.numeric(newdata[j, ])
		interp.res[j, 3] <- .cressman.pixel(xy.g, locations, values, nmin, nmax, spheric)
	}

	return(interp.res)
}

kriging.interp <- function(locations, values, newdata, vgm, nmin, nmax, spheric){
	# locations: matrix X Y
	# values: vector
	# newdata: matrix X Y
	# vgm: gstat variogram model object
	# nim: integer
	# nmax: integer
	# spheric: logical
	ngrd <- nrow(newdata)
	interp.res <- matrix(NA, nrow = ngrd, ncol = 4)
	interp.res[, 1:2] <- newdata
	for(j in seq(ngrd)){
		xy.g <- as.numeric(newdata[j, ])
		interp.res[j, 3:4] <- .okriging.pixel(xy.g, locations, values, vgm, nmin, nmax, spheric)
	}

	return(interp.res)
}

spheremap.interp <- function(locations, values, newdata, nmin, nmax, spheric){
	ngrd <- nrow(newdata)
	interp.res <- matrix(NA, nrow = ngrd, ncol = 3)
	interp.res[, 1:2] <- newdata
	for(j in seq(ngrd)){
		xy.g <- as.numeric(newdata[j, ])
		interp.res[j, 3] <- .spheremap.pixel(xy.g, locations, values, nmin, nmax, spheric)
	}

	return(interp.res)
}

##########################################################################################

.idw.pixel <- function(xy, XY, values, nmin, nmax, spheric, p){
	# xy: vector
	# values: vector
	# XY: matrix (X, Y)
	# nmin, nmax: integer
	# spheric: logical
	# p: real
	ldst <- .distance.pixel(xy, XY, nmin, nmax, spheric)
	npts <- ldst$npts
	dst <- ldst$dst[1:npts]
	z.s <- values[ldst$order][1:npts]

	Wk <- 1/dst^p
	if(any(is.infinite(Wk))){
		wi <- Wk[!is.infinite(Wk)]
		wi <- 1 - ((max(wi) - wi)/diff(range(wi)))
		Wk[!is.infinite(Wk)] <- wi/2
		Wk[is.infinite(Wk)] <- 1
	}
	sum(Wk * z.s) / sum(Wk)
}

.shepard.pixel <- function(xy, XY, values, nmin, nmax, spheric, p){
	# xy: vector
	# values: vector
	# XY: matrix (X, Y)
	# nmin, nmax: integer
	# spheric: logical
	# p: real
	ldst <- .distance.pixel(xy, XY, nmin, nmax, spheric)
	npts <- ldst$npts
	rmax <- ldst$rmax
	dst <- ldst$dst[1:npts]
	z.s <- values[ldst$order][1:npts]

	Wk <- ((rmax - dst) / (rmax * dst))^p
	if(any(is.infinite(Wk))){
		wi <- Wk[!is.infinite(Wk)]
		wi <- 1 - ((max(wi) - wi)/diff(range(wi)))
		Wk[!is.infinite(Wk)] <- wi/2
		Wk[is.infinite(Wk)] <- 1
	}
	sum(Wk * z.s) / sum(Wk)
}

.barnes.pixel <- function(xy, XY, values, nmin, nmax, spheric, p){
	# xy: vector
	# values: vector
	# XY: matrix (X, Y)
	# nmin, nmax: integer
	# spheric: logical
	# p: real (smoothing parameter) ranges between 0 and 1
	ldst <- .distance.pixel(xy, XY, nmin, nmax, spheric)
	npts <- ldst$npts
	rmax <- ldst$rmax
	dst <- ldst$dst[1:npts]
	z.s <- values[ldst$order][1:npts]

	Wk <- exp(-(dst/(p * rmax))^2)
	sum(Wk * z.s) / sum(Wk)
}

.cressman.pixel <- function(xy, XY, values, nmin, nmax, spheric){
	# xy: vector
	# values: vector
	# XY: matrix (X, Y)
	# nmin, nmax: integer
	# spheric: logical
	ldst <- .distance.pixel(xy, XY, nmin, nmax, spheric)
	npts <- ldst$npts
	rmax <- ldst$rmax
	dst <- ldst$dst[1:npts]
	z.s <- values[ldst$order][1:npts]

	Wk <- (rmax^2 - dst^2)/(rmax^2 + dst^2)
	sum(Wk * z.s) / sum(Wk)
}

.okriging.pixel <- function(xy, XY, values, vgm, nmin, nmax, spheric){
	# Ordinary Kriging
	# xy: vector
	# values: vector
	# XY: matrix (X, Y)
	# vgm: gstat variogram model object 
	# nmin, nmax: integer
	# spheric: logical
	ldst <- .distance.pixel(xy, XY, nmin, nmax, spheric)
	npts <- ldst$npts
	# dst <- ldst$dst[ldst$order][1:npts]
	z.s <- values[ldst$order][1:npts]

	xy.s <- XY[ldst$order, ][1:npts, , drop = FALSE]

	dS <- distance.matrix(xy.s, xy.s, spheric)
	np <- nrow(dS)
	dG <- distance.vector(xy, xy.s, spheric)

	vS <- vgm$psill[2] - predict.vgm(dS, vgm)
	vS <- rbind(cbind(vS, 1), c(rep(1, np), 0))
	ivS <- inverse.matrix(vS)

	vG <- vgm$psill[2] - predict.vgm(dG, vgm)
	vG <- c(vG, 1)

	W <- ivS %*% vG
	var1.pred <- sum(z.s * W[1:np, 1])
	var1.var <- vgm$psill[2] - sum(W[ ,1] * vG)
	c(var1.pred, var1.var)
}

.spheremap.pixel <- function(xy.g, XY, values, nmin, nmax, spheric){
	# xy.g: vector
	# values: vector
	# XY: matrix (X, Y)
	# nmin, nmax: integer
	# spheric: logical

	ldst <- .distance.pixel(xy.g, XY, nmin, nmax, spheric)
	npts <- ldst$npts
	rmax <- ldst$rmax
	dst <- ldst$dst[1:npts]
	dat.s <- values[ldst$order]
	z.s <- dat.s[1:npts]
	xy.s <- XY[ldst$order, ][1:npts, , drop = FALSE]

	if(all(dst < 1e-10)) return(sum(z.s)/npts)

	Sk <- rep(0, npts)
	crt1 <- dst <= rmax/3
	crt2 <- dst > rmax/3 & dst <= rmax
	Sk[crt1] <- 1/dst[crt1]
	Sk[crt2] <- ((dst[crt2]/rmax) - 1)^2 * (27/(4 * rmax))
	if(any(is.infinite(Sk))){
		sk <- Sk[!is.infinite(Sk)]
		sk <- 0.1 * (sk - min(sk))/(max(sk) - min(sk))
		Sk[!is.infinite(Sk)] <- sk
		Sk[is.infinite(Sk)] <- 1
	}

	Tk <- rep(NA, npts)
	Wk <- rep(NA, npts)
	for(k in 1:npts){
		if(spheric){
			cosd.kl <- cos.spheric.distance(xy.s[k, ], xy.s[-k, , drop = FALSE])
			cosd.jk <- cos(dst[k])
			cosd.jl <- cos(dst[-k])
			sind.jk <- sin(dst[k])
			sind.jl <- sin(dst[-k])
			cos_theta.kl <- (cosd.kl - cosd.jk * cosd.jl)/(sind.jk * sind.jl)
		}else{
			xx.kl <- (xy.s[k, 1] - xy.g[1]) * (xy.s[-k, 1] - xy.g[1])
			yy.kl <- (xy.s[k, 2] - xy.g[2]) * (xy.s[-k, 2] - xy.g[2])
			cos_theta.kl <- (xx.kl + yy.kl)/(dst[k] * dst[-k])
		}
		Tk[k] <- sum(Sk[-k] * (1 - cos_theta.kl), na.rm = TRUE)
		Wk[k] <- Sk[k]^2 * (1 + Tk[k] /sum(Sk[-k]))
	}

	if(spheric){
		dlon.jk <- (xy.g[1] - xy.s[, 1]) * cos(xy.g[2] * pi/180)
		dlat.jk <- xy.g[2] - xy.s[, 2]
	}else{
		dlon.jk <- xy.g[1] - xy.s[, 1]
		dlat.jk <- xy.g[2] - xy.s[, 2]
	}

	deltaZ.k <- 0
	if(all((z.s[1] - z.s) != 0)){
		dz.lon.k <- rep(NA, npts)
		dz.lat.k <- rep(NA, npts)
		for(k in 1:npts){
			dst.kl <- distance.vector(xy.s[k, ], xy.s[-k, ], spheric)
			if(spheric){
				dlon.lk <- (xy.s[-k, 1] - xy.s[k, 1]) * cos(xy.s[-k, 2] * pi/180)
				dlat.lk <- xy.s[-k, 2] - xy.s[k, 2]
				dst.kl <- 6378.388 * dst.kl
			}else{
				dlon.lk <- xy.s[-k, 1] - xy.s[k, 1]
				dlat.lk <- xy.s[-k, 2] - xy.s[k, 2]
			}
			dz.lon.k[k] <- sum(Wk[-k] * (z.s[-k] - z.s[k]) * dlon.lk * (1/dst.kl^2))/sum(Wk[-k]) 
			dz.lat.k[k] <- sum(Wk[-k] * (z.s[-k] - z.s[k]) * dlat.lk * (1/dst.kl^2))/sum(Wk[-k])
		}

		v <- 0.1 * diff(range(dat.s, na.rm = TRUE)) / max(sqrt(dz.lon.k^2 + dz.lat.k^2))
		deltaZ.k <- (dz.lon.k * dlon.jk + dz.lat.k * dlat.jk) * (v/(v + dst))
	}

	sum(Wk * (z.s + deltaZ.k))/sum(Wk)
}

##########################################################################################

predict.vgm <- function(x, vgm){
	# x: vector or matrix
	# vgm: gstat variogram model object 
	Nug <- vgm$psill[1]
	Sill <- vgm$psill[2]
	Rg <- vgm$range[2]
	switch(as.character(vgm$model[2]),
			"Gau" = Nug + (Sill - Nug) * (1 - exp(-3 * (x/Rg)^2)),
			"Exp" = Nug + (Sill - Nug) * (1 - exp(-3 * x/Rg)),
			"Sph" = ifelse(x <= Rg, Nug + (Sill - Nug) * (1.5 * (x/Rg) - 0.5 * (x/Rg)^3), Sill),
			"Pen" = ifelse(x <= Rg, Nug + (Sill - Nug) * (1.875 * (x/Rg) - 1.25 * (x/Rg)^3 + 0.375 * (x/Rg)^5), Sill)
		)
}

inverse.matrix <- function (X, tol = 1e-13)
{
	# Pseudoinverse
	# http://web.cs.iastate.edu/~cs577/handouts/svd.pdf
	s <- svd(X)
	e <- s$d
	e[e > tol] <- 1/e[e > tol]
	return(s$v %*% diag(e) %*% t(s$u))
}
