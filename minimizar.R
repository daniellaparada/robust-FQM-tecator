## -----------------------------------------------------------
## Classical and robust estimation under a quadratic model
## -----------------------------------------------------------
## INPUT
##
## yy      : response
## xx_coef : estimated fourier coefficients
## freq    : frequency cut-off
## fLoss   : loss function ('ls','lmrob')
## cterho  : when using 'lmrob' tuning constant for the M-step
## 	     : default 4.685061  # use 3.443689 for lower efficiency 
## nresamp : default 5000
 
##
## -----------------------------------------------------------
## OUTPUT
##
## slope  : slope parameters
## gamma  : gamma parameters
## value  : minimization value
## scale  : estimated scale (if used)
## conv   : convergence flag
## -----------------------------------------------------------

minimizar <- function (yy, xx_coef, zeta, freq, fLoss,  cterho=4.685061, nresamp=5000) {

    ## Initialize
    cv   <- NA
    tr   <- NA
    aicM <- NA
    
    ## Design matrix
    X <-  cbind(xx_coef, zeta)
    
    ## Minimization menu
    if (fLoss == 'ls') {
        fit <- lm(yy ~ X)
        cf  <- fit$coef #includes intercept in position  1
        vv  <- 1  
        ss  <- sqrt(mean(fit$res^2))
    }
  
    if (fLoss == 'lmrob') {
        control <- lmrob.control(trace.level = 0,         #  
                                 nResample   =  nresamp,  #  
                                 tuning.psi = cterho,     #  
                                 subsampling = 'simple',  #
                                 rel.tol     = 1e-5,      #  
                                 refine.tol  = 1e-5,      # 
                                 k.max       = 2e3,       #  
                                 maxit.scale = 2e3,       #  
                                 max.it      = 2e3)       #  
        fit  <- lmrob(yy ~ X, control = control)
        if (fit$init.S$converged) {
            cf <- fit$coef  #includes intercept in position  1
            ss <- fit$scale
            vv <- sum(Mpsi(fit$res / ss,
                            cc  = control$tuning.psi,
                            psi = control$psi, deriv = -1))
            cv <- fit$converged
            
        } else {
            stop('No convergence for the  S-estimador. The subsample is discarded.')
        }
    }
 
    #############################
    ## Estimated parameters
    #############################

    aa <- cf[-1] 
    slope_par <- aa[1:freq]
    ordenada <- cf[1]
    gamma_par <- aa[-(1:freq)]
    
    return(list(ordenada = ordenada, slope = slope_par, gamma = gamma_par,
                value = vv, scale = ss, conv = cv))

}



## -----------------------------------------------------------
## Estimation under a linear model
## -----------------------------------------------------------
## INPUT
##
## yy      : response
## xx_coef : estimated fourier coefficients
## freq    : frequency cut-off
## fLoss   : loss function ('ls','lmrob')
##
 
##
## -----------------------------------------------------------
## OUTPUT
##
## slope  : slope parameters 
## value  : minimization value
## scale  : estimated scale (if used)
## conv   : convergence flag
## -----------------------------------------------------------

minimizar.lineal <- function (yy, xx_coef,  freq, fLoss,  cterho=4.685061, nresamp=5000) {

    ## Initialize
    cv   <- NA
    tr   <- NA
    aicM <- NA
    
    ## Design matrix
    X <-  xx_coef 
    
    ## Minimization menu
    if (fLoss == 'ls') {
        fit <- lm(yy ~ X)
        cf  <- fit$coef #includes intercept in position  1
        vv  <- 1  
        ss  <- sqrt(mean(fit$res^2))
    }
  
    if (fLoss == 'lmrob') {
        control <- lmrob.control(trace.level = 0,         # 0
                                 nResample   =  nresamp,  # 5000 default
                                 tuning.psi = cterho,     # for 85% eff use 3.443689
                                 subsampling = 'simple',  #
                                 rel.tol     = 1e-5,      #  
                                 refine.tol  = 1e-5,      #  
                                 k.max       = 2e3,       #  
                                 maxit.scale = 2e3,       #  
                                 max.it      = 2e3)       #  
        fit  <- lmrob(yy ~ X, control = control)
        if (fit$init.S$converged) {
            cf <- fit$coef  #includes intercept in position  1
            ss <- fit$scale
            vv <- sum(Mpsi(fit$res / ss,
                            cc  = control$tuning.psi,
                            psi = control$psi, deriv = -1))
            cv <- fit$converged
            
        } else {
            stop('No convergence for the  S-estimador. The subsample is discarded.')
        }
    }

    #############################
    ## Estimated parameters
    #############################
    aa <- cf[-1] 
    slope_par <- aa[1:freq]
    ordenada <- cf[1]
     
    return(list(ordenada = ordenada, slope = slope_par,  
                value = vv, scale = ss, conv = cv))
}


