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


