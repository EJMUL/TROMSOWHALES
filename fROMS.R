readROMS <- function(roms.dir='R:/Jofrid/Kaldfjord_160m_2014', date='20141101', var='u', 
                     xlim=c(18.3, 18.8), ylim=c(69.6, 69.9), to.xyz=F) {

  require(ncdf4)
  
  varname <- var

  filelist <- dir(roms.dir)
  filedates <- unlist(lapply(filelist, function(x) {
    substr(tail(unlist(strsplit(x, '_')), 1), 1, 8)
  }))  
  
  which.file <- match(date, filedates)
  nc <- nc_open(paste(roms.dir, filelist[which.file], sep='/'))
  
  dim.nam <- paste(c('lon', 'lat'), var, sep='_')
  
  if(var=='temp' | var=='salt') {
    dim.nam <- paste(c('lon', 'lat'), 'rho', sep='_')
  }
  
  eval(parse(text=paste('lon <- ncvar_get(nc, nc$var$', dim.nam[1], ')', sep='')))
  eval(parse(text=paste('lat <- ncvar_get(nc, nc$var$', dim.nam[2], ')', sep='')))
  eval(parse(text=paste('var <- ncvar_get(nc, nc$var$', var, ')', sep='')))
  
  nc_close(nc)
  
  if(any(!is.na(xlim))) {
    lon.clip <- which(lon>=xlim[1] & lon<=xlim[2])  
    lat.clip <- which(lat>=ylim[1] & lat<=ylim[2])  
  }
  
  
  if(to.xyz) {
    xyz <- data.frame(lon=lon[intersect(lon.clip, lat.clip)], 
                    lat=lat[intersect(lon.clip, lat.clip)], 
                    var=var[,,1][intersect(lon.clip, lat.clip)])
    names(xyz)[3] <- varname
    xyz
  } else {
    out.l <- list(lon=lon, lat=lat, var=var)
    names(out.l)[3] <- varname
    out.l
  }  
}

cropROMS <- function(uc=u, vc=v, xdim=c(100:200), ydim=c(1:150)) {
  vdim=list(xdim=c(xdim, max(xdim)+1), ydim=ydim[-length(ydim)])

  uc$lon <- uc$lon[xdim, ydim]
  uc$lat <- uc$lat[xdim, ydim]
  uc$u <- uc$u[xdim, ydim,]

  vc$lon <- vc$lon[vdim$xdim, vdim$ydim]
  vc$lat <- vc$lat[vdim$xdim, vdim$ydim]
  vc$v <- vc$v[vdim$xdim, vdim$ydim,]
  list(u=uc, v=vc)
}


imageROMS <- function(u=u, v=v, arr.dens=2, scale=100, arr.col='slategrey') {
  require(fields)
  options(warn=-1)
  del.dims <- c(max(c(dim(u$u)[1], dim(v$v)[1])),
               max(c(dim(u$u)[2], dim(v$v)[2])))
  im.dims <- list(c(1:(del.dims[1]-1)), 
                  c(1:(del.dims[2]-1)))
  
  arr.starts <- as.matrix(expand.grid(seq(1, del.dims[1]-1, by=arr.dens),
                   seq(1, del.dims[2]-1, by=arr.dens), KEEP.OUT.ATTRS = F))
  arr.stops <- t(apply(arr.starts, 1, function(x) {
    c(x[1]+(scale*u$u[x[1],x[2],1]), x[2]+(scale*v$v[x[1], x[2], 1]))
  }))
  image.plot(im.dims[[1]], 
             im.dims[[2]], 
             sqrt(u$u[,-del.dims[2],1]^2+v$v[-del.dims[1],,1]^2), 
             col=colorRampPalette(brewer.pal(3, 'Reds'))(256),
             axes=F, xlab='', ylab='')
              box()
  arrows(arr.starts[,1], arr.starts[,2], arr.stops[,1], arr.stops[,2], length=0.05, col=arr.col)
  options(warn=0)
}


## Example: 
dates <- as.character(c(20141101:20141130))
for(i in 1:length(dates)) {
  cat('Reading Norfjords for', dates[i], '\n')
  flush.console()
  u <- readROMS(date=dates[i])
  v <- readROMS(date=dates[i], var='v')
  uv <- cropROMS()
  png(paste('Norfjords', dates[i], 'png', sep='.'), width=20, height=20, units='cm', res=500, pointsize=10)
  imageROMS(uv$u, uv$v)
  mtext(side=3, line=1, cex=1.3, format(strptime(dates[i], '%Y%m%d'), '%d-%b-%Y'))
  dev.off()
}

