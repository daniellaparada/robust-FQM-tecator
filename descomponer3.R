## -------------------------------------------------------
## Computes an estimate for the covariance operator
## and its eigenfunctions
## -------------------------------------------------------
##
## xx     : matrix n by k (n is the sample size), it has to be centered first
## dt     : mesh size
## method : 'cl'   classical
##          'ger'  robust (Gervini, 2008)
##          'gerS' robust (Gervini, 2008) with S-ordering 
## -------------------------------------------------------


## Center the data
centrar <- function(xx, method) {
    nn <- nrow(xx)
    cc <- switch(method,
                 cl = colMeans(xx),
                 pcaPP::l1median(xx))
    xcenter <- t(t(xx) - cc)
    return(list(centro = cc, xcenter = xcenter))     
}

###################################################
# COMPUTES THE EIGENVALUES AND EIGENFUNCIONS 
###################################################

descomponer3 <- function (xx, dt, method) {
	#############################
	# Center the data
	#############################

	salida=centrar(xx,method)
	xcenter=salida$xcenter
	#############################################
	# Computes the eigenfunctions and eigenvalues
	##############################################

	switch(method,
           cl   = pc_cl2(xcenter, dt),
           ger  = pc_esf2(xcenter, dt, auto = FALSE),
           gerS = pc_esf2(xcenter, dt, auto = TRUE))
}
 
###################################################
## Classical principal components and eigenvalues
###################################################

pc_cl2 <- function(xx, dt) {
    ## Eigen
    cov_mat <- crossprod(xx) / nrow(xx)
    cov_dec <- eigen(cov_mat, symmetric = TRUE)
    pc_cl   <- cov_dec$vectors
    norms   <- sqrt(colSums(pc_cl^2) * dt)
    pc_cl_n <- pc_cl %*% diag(1 / norms)
    return(list(autofun=pc_cl_n, autoval=cov_dec$values))
}

###################################################
## Spherical principal components and eigenvalues
###################################################

pc_esf2 <- function(xx, dt,  auto) {
    ## Spherical data
    ww <- sqrt(rowSums(xx^2))
    ww <- replace(ww, which(ww <= 1e-50), 1e-50)
    xe <- xx / ww
    ## Eigen
    nn       <- nrow(xx)
    mm       <- colMeans(xe)
    cov_esf  <- crossprod(xe) / nn
    cov_dec  <- eigen(cov_esf, symmetric = TRUE) 
    pc_esf   <- cov_dec$vectors
    av_esf   <- cov_dec$values
    norms    <- sqrt(colSums(pc_esf^2) * dt)
    pc_esf_n <- pc_esf %*% diag(1 / norms)
    ## Sort by spherical eigenvalues 
    if (auto) {
        av_esf   <- apply(xx %*% pc_esf_n * dt, 2, s_scale)^2
        ordendec <- order(av_esf, decreasing = TRUE)
        pc_esf_n <- pc_esf_n[, ordendec]
	av_esf  <- av_esf[ordendec]
    }
    return(list(autofun=pc_esf_n, autoval=av_esf))
}

##############################
## Version of FastS
##############################

s_scale <- function(u, b = 0.5 , cc = 1.54764) {
    ## find the scale, full iterations
    max.it <- 200
    sc     <- median(abs(u)) / 0.6745
    i      <- 0
    eps    <- 1e-20
    ## magic number alert
    err <- 1
    while(((i <- i + 1) < max.it) && (err > eps)) {
        sc2 <- sqrt(sc^2 * mean(Mchi(u / sc, cc, psi = 'bisquare')) / b)
        err <- abs(sc2 / sc - 1)
        sc  <- sc2
    }
    return(sc)
}
