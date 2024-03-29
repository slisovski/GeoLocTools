setupGeolocation <- function(force = NULL) {

  reqPackages <- c("devtools","digest","geosphere","raster","fields","forecast", "mapdata",
                   "circular","truncnorm","parallel","bit","rgdal","CircStats","Rcpp", "grid",
                   "RcppArmadillo","ggmap","ggsn","sp","maptools","rgeos","MASS", "svglite", "terra")

  ### CRAN
  get.packages <- reqPackages[!(reqPackages%in%installed.packages()[,"Package"])]

  if(length(get.packages)>0){
    install.packages(get.packages, repos = "https://cloud.r-project.org/", dependencies = TRUE, force = T)
  }
  data(wrld_simpl)

  null <- suppressMessages(lapply(reqPackages, require, character.only = TRUE))

  # Install necessary packages from Github using the devtools library #

  reqGitHPackages <- c("SGAT", "TwGeos", "PolarGeolocation", "GeoLight", "FLightR", "probGLS")
  get.packages    <- reqGitHPackages[!(reqGitHPackages %in% installed.packages()[,"Package"])]

  if(any("GeoLight"%in%get.packages)) {
    cat("Install and load R package GeoLight from GitHub.")
    suppressMessages({
    install_github("SLisovski/GeoLight", quiet = T)
    suppressMessages(library(GeoLight))})
    cat(" ..... Done!\n")
  } else suppressMessages(library(GeoLight))
  if(any("SGAT"%in%get.packages)) {
    cat("Install and load R package SGAT from GitHub.")
    suppressMessages({
      install_github("SWotherspoon/SGAT", quiet = T)
      suppressMessages(library(SGAT))})
    cat(" ..... Done!\n")
  } else suppressMessages(library(SGAT))
  if(any("probGLS"%in%get.packages)) {
    cat("Install and load R package probGLS from GitHub.")
    suppressMessages({
      install_github("benjamin-merkel/probGLS", quiet = T)
      suppressMessages(library(probGLS))})
    cat(" ..... Done!\n")
  } else suppressMessages(library(probGLS))
  if(any("TwGeos"%in%get.packages)) {
    cat("Install and load R package TwGeos from GitHub.")
    suppressMessages({
      install_github("SLisovski/TwGeos", quiet = T)
      suppressMessages(library(TwGeos))})
    cat(" ..... Done!\n")
  } else suppressMessages(library(TwGeos))
  if(any("PolarGeolocation"%in%get.packages)) {
    cat("Install and load R package PolarGeolocation from GitHub.")
    suppressMessages({
      install_github("SLisovski/PolarGeolocation", quiet = T)
      suppressMessages(library(PolarGeolocation))})
    cat(" ..... Done!\n")
  } else suppressMessages(library(PolarGeolocation))
  if(any("FLightR"%in%get.packages)) {
    cat("Install and load R package FLightR from GitHub.")
    suppressMessages({
      install_github("eldarrak/FLightR", quiet = T)
      suppressMessages(library(FLightR))})
    cat(" ..... Done!\n")
  } else suppressMessages(library(FLightR))
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

world <- rnaturalearth::ne_countries(
  scale = "small",
  returnclass = "sf")

if(pacific) {
  world <- world %>%
    st_break_antimeridian(lon_0 = 180) %>% # insert this before transformation
    st_transform(crs = "+proj=longlat +lon_0=180")
}

nrows <- abs(lat[2L] - lat[1L]) / cell.size
ncols <- abs(lon[2L] - lon[1L]) / cell.size

grid <- rast(
  nrows = nrows,
  ncols = ncols,
  xmin = min(lon),
  xmax = max(lon),
  ymin = min(lat),
  ymax = max(lat),
  crs = crs(world)
)

grid <- rasterize(world, grid, field = 1, silent = TRUE)
grid <- is.na(grid)

switch(mask,
       sea = {},
       land = {
         grid <- subst(grid, from = c(0,1), to = c(1,0))
       },
       none = {
         grid <- subst(grid, from = c(0,1), to = c(1,1))
       }
)
return(grid)
}


landMask <- function(xlim, ylim, cell.size = 1, land = TRUE, pacific = FALSE) {
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

maskStack <- function(xlim, ylim, behav, n = 2) {
  # This ensures that the bird can fly over sea during the movement period,
  # but that when it is stationary, it must be on land. To do so, we create a
  # raster where the bird can go anywhere, and one where it is restricted to land.
  # Then we tell it to lookup any given lat lon on either depending on whether
  # the bird is moving or not.

  #create an empty raster
  r <- raster(nrows = n * diff(ylim), ncols = n * diff(xlim), xmn = xlim[1],
              xmx = xlim[2], ymn = ylim[1], ymx = ylim[2], crs = proj4string(wrld_simpl))
  ## land for stationary periods
  rs <- cover(rasterize(elide(wrld_simpl, shift = c(-360, 0)), r, 1, silent = TRUE),
              rasterize(wrld_simpl, r, 1, silent = TRUE),
              rasterize(elide(wrld_simpl,shift = c(360, 0)), r, 1, silent = TRUE))
  # anywhere for migratory periods
  rm <- rs
  rm[] <- 1

  st    <- stack(rs, rm)
  index <- ifelse(behav, 1, 2)

  list(index = index, mask = st)
}



maskExtractor <- function(s) {

  xbin  <- seq(xmin(s$mask),xmax(s$mask),length=ncol(s$mask)+1)
  ybin  <- seq(ymin(s$mask),ymax(s$mask),length=nrow(s$mask)+1)
  mask  <- as.array(s$mask)[nrow(s$mask):1,,sort(unique(s$index)),drop=FALSE]
  index <- s$index

  function(p) mask[cbind(.bincode(p[,2],ybin),.bincode(p[,1],xbin),index)]
}

