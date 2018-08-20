setupGeolocation <- function(forceGeoLight = FALSE) {

  reqPackages <- c("devtools","digest","geosphere","raster","fields","forecast",
                   "circular","truncnorm","parallel","bit","rgdal","CircStats","Rcpp", "grid",
                   "RcppArmadillo","ggmap","ggsn","sp","maptools","rgeos","MASS")

  ### CRAN
  get.packages <- reqPackages[!(reqPackages %in% installed.packages()[,"Package"])]

  if(length(get.packages)>0){
    install.packages(get.packages,repos = "https://cloud.r-project.org/", dependencies = TRUE)
  }

  library(devtools)
  library(GeoLight)
  library(maptools); data(wrld_simpl)
  library(MASS)

  message("All nessesary packges available from CRAN are installed!")


  # Install necessary packages from Github using the devtools library #

  reqGitHPackages <- c("SGAT","TwGeos", "PolarGeolocation","GeoLight", "FLightR")
  get.packages    <- reqGitHPackages[!(reqGitHPackages %in% installed.packages()[,"Package"])]

  library(devtools)
  if(any("SGAT"%in%get.packages)) install_github("SWotherspoon/SGAT"); library(SGAT)
  if(any("TwGeos"%in%get.packages)) install_github("SLisovski/TwGeos") ; library(TwGeos)
  if(any("PolarGeolocation"%in%get.packages)) install_github("SLisovski/PolarGeolocation"); library(PolarGeolocation)
  if(forceGeoLight || any("GeoLight"%in%get.packages)) install_github("SLisovski/GeoLight", ref = "Update_2.01", force = T); library(GeoLight)
  if(any("FLightR"%in%get.packages)) install_github("eldarrak/FLightR"); library(FLightR)

  message("You are all set!")

}


findHEZenith <- function(twl, tol = 0.08, range=c(250,400)){

  z <- seq(89, 99, by = 0.25)

  lats  <- apply(cbind(z), 1, function(x) thresholdPath(twl$Twilight, twl$Rise, zenith = x, tol=tol)$x[,2])
  sds   <- apply(lats[range[1]:range[2],], 2, sd, na.rm = T)

  colsT <- data.frame(sd = seq(min(sds)-0.1, max(sds)+0.1, length = 100), col = heat.colors(100))

  opar <- par(mfrow = c(2,1), mar=c(4,4,1,1))
  matplot(lats, col = as.character(colsT$col[cut(sds, breaks = colsT[,1], labels = F)]),
          lty = 1, type = "l", xlab = "", ylab = "Latitude", las = 1, xaxt = "n")
  lines(lats[,which.min(sds)], lwd = 3)
  abline(v = range, lty = 2, col = "cornflowerblue")
  axis(1, at = seq(1, nrow(twl), length = 5),
       labels = format(as.POSIXct(seq(twl$Twilight[1], twl$Twilight[nrow(twl)], length = 5)), "%d-%b"))
  plot(z, sds, las = 1, type = "o", pch = 16, cex = 1.3, col = as.character(colsT$col[cut(sds, breaks = colsT[,1], labels = F)]),
       xlab = "zenith", ylab = "sd in latitude (within range)")

  abline(v = z[which.min(sds)], lty = 3)
  mtext(paste("zenith =", z[which.min(sds)]), 3, line = -1.5, cex = 1.1)
  par(opar)

  return(z[which.min(sds)])

}



makeGroups <- function(grouped){
  g <- rep(0, length(grouped))
  while(any(g<1)) {
    ind1 <- min(which(g<1))
    if(length(unique(grouped[ind1:length(g)]))==1) {
      ind2 <- length(g)
    } else {
      ind2 <- min(which(grouped[ind1:length(g)]!=grouped[ind1])-1)+(ind1-1)
    }
    if(grouped[ind1]) {
      g[ind1:ind2] <- rep(max(g, na.rm = T)+1, length(ind1:ind2))
    } else {
      g[ind1:ind2] <- max(g, na.rm = T)+seq_len(length(ind1:ind2))
    }
  }
  return(g)
}



makeGrid <- function(lon = c(-180, 180), lat = c(-90, 90), cell.size = 1, mask = "sea", pacific = FALSE) {
  data(wrld_simpl, package = "maptools", envir = environment())
  if(pacific){
    wrld_simpl <- nowrapRecenter(wrld_simpl, avoidGEOS = TRUE)}
  nrows <- abs(lat[2L] - lat[1L]) / cell.size
  ncols <- abs(lon[2L] - lon[1L]) / cell.size
  grid <- raster(
    nrows = nrows,
    ncols = ncols,
    xmn = min(lon),
    xmx = max(lon),
    ymn = min(lat),
    ymx = max(lat),
    crs = proj4string(wrld_simpl)
  )
  grid <- rasterize(wrld_simpl, grid, 1, silent = TRUE)
  grid <- is.na(grid)
  switch(mask,
         sea = {},
         land = {
           grid <- subs(grid, data.frame(c(0,1), c(1,0)))},
         none = {
           grid <- subs(grid, data.frame(c(0,1), c(1,1)))
         }
  )
  return(grid)
}


land.mask <- function(xlim, ylim, cell.size = 1, land = TRUE, pacific = FALSE) {
  r <- makeGrid(xlim, ylim, cell.size, ifelse(land, "land", "sea"), pacific = pacific)
  r <- as.matrix(is.na(r))[nrow(r):1, ]
  if (land)
    r <- !r
  xbin <- seq(xlim[1], xlim[2], length = ncol(r) + 1)
  ybin <- seq(ylim[1], ylim[2], length = nrow(r) + 1)

  function(p) {
    r[cbind(.bincode(p[, 2], ybin), .bincode(p[, 1], xbin))]
  }
}

