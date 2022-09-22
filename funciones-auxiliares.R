## -------------------------------------------------------
## 
## -------------------------------------------------------
## INPUT
##
## y	       :response
## xcenter   : covariates
## ttt       : grid (equiespaced)
## freq      : Number of eigenfunctions to be considered
## varexp    : Percentage of variance explained, if varexp=1 the value of freq is taken as the number of PC

## cov_type  : 'cl' classic 
##             'gerS' robust (Spherical PC)
 
## fLoss     : loss function 
## ajuste    : 'lineal' fits a linear regression model
##           : 'quadra' fits a quadratic one


 
  
estimar <- function (y,xcenter,  ttt, ajuste='lineal', freq=NULL,  cov_type,  fLoss,   cterho=3.443689, nresamp=5000, varexp=varexp){
  
  ## Some samples may be discarded...
   nn=length(y)
   lt=length(ttt)
   dt=ttt[2]-ttt[1]
  
     

  ################################
  # Decomposition and projection #
  #####################
  cov_dec1 <- descomponer3(xcenter, dt, cov_type)
  autofun <- cov_dec1$autofun
  autoval <- cov_dec1$autoval
  	 
  
  if(varexp!=1){
  	porcen_acum <- cumsum(autoval)/sum(autoval)

  	##############################################
  	# Which are >= fixed % of the total variation
  	##############################################

  	cuales <- 1*(varexp <= porcen_acum)
  	names(cuales) <- 1: length(autoval)
  	distintoscero <- as.numeric(names(cuales[cuales!=0]))
  
  	freq= min(distintoscero)
  	}

  porcentaje= sum(autoval[1:freq])/sum(autoval)

  cov_dec=matrix(NA, ncol=freq,nrow=lt)
  cov_dec[, 1:freq] <- autofun[, 1:freq] 

  ## Estimated coefficients (by row)
  xx_coef <- xcenter %*% cov_dec * dt
  
  if(ajuste=='lineal'){
	est <-  minimizar.lineal(y, xx_coef,   freq, fLoss, cterho=cterho, nresamp=nresamp) 
 	######################################
      ## Slope function estimate
	######################################
     	
	est_slope_fun <- cov_dec %*% est$slope
	
	######################################
	## Intercept
	######################################

 	mu <- est$ordenada

 	######################################
	## NO COEFFICIENTS FOR THE QUADRATIC TERM
	######################################

      gamma_coef <- matrix(NA,freq,freq)
	est_gamma_fun<- NA
	}


  if(ajuste=='quadra'){
  	######################################
	## Estimated coefficients of the quadratic term
	######################################

  	dimensionz <- freq*(freq+1)/2
  	zeta <- matrix(NA, nn, dimensionz)
  	for (i in 1:nn){
    		x_aux <- xx_coef[i,]
    		z_aux <- xx_coef[i,] %*% t(xx_coef[i,])
    		zeta[i,] <- upperTriangle(z_aux, diag=TRUE, byrow=TRUE)
  		}
    
  	#########################
  	## Parameter estimation #
  	#########################
  	est <-  minimizar(y, xx_coef, zeta, freq, fLoss, cterho=cterho, nresamp=nresamp) 
 
	######################################
      ## Slope function estimate
	######################################
 
 	est_slope_fun <- cov_dec %*% est$slope
      
 	######################################
  	## Intercept
	######################################

 	mu <- est$ordenada
       
 	######################################
 	## Start Quadratic coefficients estimates
	######################################

  	gamma_coef <- matrix(NA,freq,freq)

	###############################################################
	# Fill the upper triangle and diagonal with the estimates 
	##############################################################

  	upperTriangle(gamma_coef,diag=TRUE,byrow=TRUE) <- est$gamma 
  	
	##################################################
	# Split the upper triangle (without the diagonal) and 
	# divide by  2 to obtain the estimated coefficients
	##################################################

	aux <- upperTriangle(gamma_coef, diag=FALSE,byrow=TRUE)/2 
  	
	##################################################
	#Redefine the upper triangle (without the diagonal) 
	# with the obtained coefficients 
	##################################################

	upperTriangle(gamma_coef, diag=FALSE,byrow=TRUE) <- aux 
  	
	##################################################
	# Fill the lower triangle ((without the diagonal) 
	# with these coefficients 
	##################################################

	lowerTriangle(gamma_coef, diag=FALSE,byrow=FALSE) <- aux 

  	#####################################################################
	# Obtain the estimated kernel of the quadratic operator as
  	# sum_{i} gamma[i,i] phi_i(t) phi_i (s)  + 2 sum_{i\ne j} gamma[i,j] phi_i(t) phi_j (s)
  	# over the grid points  which are stored in cov_dec[,i] and cov_dec[,j]
  	# leading to a  lt*lt matrix
	# to ensure symmetry compute (vhat(t,s)+vhat(s,t))/2
  	#####################################################################
 
  	est_gamma_fun <- cov_dec %*% gamma_coef %*% t(cov_dec) 
  	est_gamma_fun <- (est_gamma_fun +t(est_gamma_fun ) )/2

		}

  	#################
  	## Save results #
  	#################	
  	return(list(alfa=mu, beta=est_slope_fun , gamma=est_gamma_fun, 
                  converge=est$conv, sigma=est$scale,freq=freq, porcentaje=porcentaje,
                  autofun=cov_dec, gama_coef=gamma_coef, slope_coef=est$slope))
       
  	}
     

## ---------------------------------------------------------------------
## Functions for the real data set
##  
## ---------------------------------------------------------------------


submuestra <- function (cuales, covariable) {
    ##################################################
    ## Construyo la submuestra (para train/test)
    ##################################################

    y  <- tecator$y$Fat[cuales]	       
    ab <- tecator$absorp.fdata
    t  <- ab$argvals
    nn <- length(y)
      
    if (covariable == 'd0') {
        x        <- ab$data[cuales, ]
    }
  
    if (covariable == 'd1') {
        ab1      <- fdata.deriv(ab, nderiv = 1) # primera derivada
        ab1$data <- ab1$data[cuales, ]
        x        <- ab1$data
    }
    if (covariable == 'd2'){
        ab2      <- fdata.deriv(ab, nderiv = 2) # segunda derivada
        ab2$data <- ab2$data[cuales, ]
        x        <- ab2$data
    }
    return(list(y = y, x = x,   t = t,   n = nn))
}
  

