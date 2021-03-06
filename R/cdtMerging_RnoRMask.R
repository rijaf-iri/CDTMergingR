
rain_no_rain.mask <- function(locations.stn, newgrid, pars.RnoR)
{
	wet.day <- pars.RnoR$wet.day
	if(wet.day <= 0) wet.day <- wet.day + 1e-13
	rnr.grd <- ifelse(newgrid@data$grd < wet.day, 0, 1)
	ij <- over(locations.stn, as(newgrid, "SpatialPixels"))
	locations.stn$rnr.stn <- ifelse(locations.stn$stn < wet.day, 0, 1)
	locations.stn$rnr.grd <- ifelse(rnr.grd[ij] < wet.day, 0, 1)
	locations.stn <- locations.stn[!is.na(locations.stn$rnr.grd), ]

	glm.binom <- glm(rnr.stn ~ rnr.grd, data = locations.stn, family = binomial(link = "logit"))

	nlon <- newgrid@grid@cells.dim[1]
	nlat <- newgrid@grid@cells.dim[2]
	rnr <- matrix(1, ncol = nlat, nrow = nlon)
	if(!is.na(glm.binom$coef[2])){
		locations.stn$rnr.res <- residuals(glm.binom)
		rnr.trend <- predict(glm.binom, newdata = newgrid, type = 'link')

		rnr.res.grd <- gstat::krige(rnr.res~1, locations = locations.stn, newdata = newgrid, debug.level = 0)
		rnr <- rnr.trend + rnr.res.grd$var1.pred

		rnr <- exp(rnr) / (1 + exp(rnr))
		### decision boundary 0.5
		rnr[rnr >= 0.5] <- 1
		rnr[rnr < 0.5] <- 0
		rnr <- matrix(rnr, ncol = nlat, nrow = nlon)
		rnr[is.na(rnr)] <- 1
		if(pars.RnoR$smooth) rnr <- smooth.matrix(rnr, 2)
	}

	return(rnr)
}
