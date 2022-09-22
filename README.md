The following real data example of the implementation of robust
estimators for functional quadratic regression models. These robust
estimators of the principal directions with robust regression estimators
based on a bounded loss function and a preliminary residual scale
estimator. This is part of a work in progress done in collaboration with
Prof. Dr. Graciela Boente.

Let’s first load some <code>R</code> packages.

    library('fda')          
    library('robustbase')   # to use lmrob 
    library(fda.usc)        # tecator dataset and fda tools
    library(gdata)          # to use uppertriangle
    library(lattice)        # to plot

Let’s now load some custom <code>R</code> functions.

    source('funciones-auxiliares.R') 
    source('descomponer3.R')         # for covariance decomposition
    source('minimizar.R')            # minimization menu

We will use the Tecator data set available in <code>R</code>
(<http://lib.stat.cmu.edu/datasets/tecator>). Each observation consists
of a spectrometric curve that corresponds to the absorbance measured on
an equally spaced grid of 100 wavelengths between 850 and 1050 nm. The
contents of fat protein and moisture were also recorded through analytic
chemistry methods. The goal of the analysis is to predict the fat
content (*y*) using some characteristics of the spectrometric curve. In
particular, we consider the first derivative (*X*) of the spectrometric
curve.

    datos <- data(tecator)
    absorp <- tecator$absorp.fdata
    absorp1 <- fdata.deriv(absorp,nderiv = 1)  #  computes the first derivative

Plots of the spectrometric data and it’s first derivative are displayed
below.

![](README_files/figure-markdown_strict/initial%20plots-1.png)![](README_files/figure-markdown_strict/initial%20plots-2.png)![](README_files/figure-markdown_strict/initial%20plots-3.png)

Functional boxplots of the spectrometric data and it’s first derivative
are displayed below.

    # First derivative
    fd1 <- fbplot(t(absorp1$data), x = absorp1$argvals, method =  "MBD",
           xlab = "Wavelength", ylab = "d(Absorbance, 1)",
           xlim = c(850, 1050), ylim = c(-0.02, 0.05))

![](README_files/figure-markdown_strict/functional%20boxplots-1.png)

    # Absorbances
    fd0<-fbplot(t(absorp$data), x = absorp$argvals, method = "MBD",
           xlab = "Wavelength", ylab = "Absorbance",
           xlim = c(850, 1050), ylim = c(2, 5.5))

![](README_files/figure-markdown_strict/functional%20boxplots-2.png)

# Robust Estimator under a Functional Linear Model

*y* = *α*<sub>0</sub> + ⟨*X*, *β*<sub>0</sub>⟩ + *ϵ*

    set.seed(124) 
    covariable <- 'd1' # first derivative  
     
    indices <- 1:215  

    # indices
    tecdatos <- submuestra(indices, covariable = covariable)
     
    indices_tecdatos <- indices
    dt_tecdatos <- tecdatos$t[2]-tecdatos$t[1]
     


    ###############################################
    ## Choose the number of principal directions
    ## equal to 4, then varexp should be 1
    ###############################################


    freq=4
    varexp=1

    est_rob_LINEAL <- estimar(y=tecdatos$y,xcenter= tecdatos$x,  
                ttt= tecdatos$t, ajuste='lineal', freq=freq,  cov_type='gerS',  
                fLoss='lmrob',cterho=3.443689, nresamp=5000, varexp=varexp) 

    #######################################
    # PERCENTAGE EXPLAINED By 4 directions
    #######################################

    est_rob_LINEAL$porcentaje

    ## [1] 0.9841754

    #0.9841754

    #############################################
    # Store the estimates and compute predictions
    # for the residuals
    #############################################

    beta_rob_LINEAL <- est_rob_LINEAL$beta
     
    phies_rob_LINEAL <- est_rob_LINEAL$autofun

    coef_rob_LINEAL <-   tecdatos$x %*% phies_rob_LINEAL * dt_tecdatos

    predichos_ROB_LINEAL <- est_rob_LINEAL$alfa + coef_rob_LINEAL%*% est_rob_LINEAL$slope_coef  


    #####################################################################################
    ## ROBUST ESTIMATOR UNDER A FQM 
    ## fat= alfa_0 + < absorp1, beta_0 > + <absorp1, Upsilon_0 absorp1>+  epsilon
    #####################################################################################

    est_rob_CUADRA <- estimar(y= tecdatos$y,xcenter= tecdatos$x,  ttt= tecdatos$t, 
                ajuste='quadra', freq=freq,  cov_type='gerS',  fLoss='lmrob',
                cterho=3.443689, nresamp=5000, varexp=varexp) 


    #######################################
    # PERCENTAGE EXPLAINED By 4 directions
    #######################################

    est_rob_CUADRA$porcentaje

    ## [1] 0.9841754

    #0.9841754


    #############################################
    # Store the estimates and compute predictions
    # for the residuals
    #############################################

    beta_rob_CUADRA <- est_rob_CUADRA$beta
    gamma_rob_CUADRA <- est_rob_CUADRA$gamma
    phies_rob <- est_rob_CUADRA$autofun

    coef_rob_CUADRA <-   tecdatos$x %*% phies_rob * dt_tecdatos

    predichos_ROB_CUADRA <- est_rob_CUADRA$alfa + coef_rob_CUADRA%*% est_rob_CUADRA$slope_coef + diag(coef_rob_CUADRA%*% est_rob_CUADRA$gama_coef%*% t(coef_rob_CUADRA))


    #########################################################
    ## Residuals from a robust fit
    #######################################################
     
    residuos_ROB_CUADRA_tecdatos <- tecdatos$y - predichos_ROB_CUADRA

    residuos_ROB_LINEAL_tecdatos <- tecdatos$y - predichos_ROB_LINEAL

    ##########################################
    # Residuals boxplots and atypical data
    #########################################

    names(residuos_ROB_LINEAL_tecdatos) <- 1:length(residuos_ROB_LINEAL_tecdatos)
    boxplot(residuos_ROB_LINEAL_tecdatos)$out

    ##          7         34         35         43         44         45        118 
    ##   5.160844 -11.582193 -17.374140 -20.565905 -18.709516 -14.669572  -7.216299 
    ##        119        121        122        126        127        128        129 
    ##  -6.140676  -6.297855  -6.917407  -5.255773 -13.431340  -7.344434 -16.054206 
    ##        130        140        143        168        171        172        185 
    ##   5.035580  -7.300928  -6.813106  -6.165244 -12.689028 -17.845070  -7.195069 
    ##        186        215 
    ## -16.402547 -16.402547

    atipicos_LINEAL_tecdatos <- as.numeric(names(boxplot(residuos_ROB_LINEAL_tecdatos)$out))

![](README_files/figure-markdown_strict/robust%20est%20FLM-1.png)

    names(residuos_ROB_CUADRA_tecdatos) <- 1:length(residuos_ROB_CUADRA_tecdatos)
    boxplot(residuos_ROB_CUADRA_tecdatos)$out

    ##         4         8        10        20        31        34        35        38 
    ## -3.385036 -1.978527  2.334157 -2.266777  2.798450 26.089185 25.532687 -2.453734 
    ##        40        43        44        45       102       108       117       121 
    ## -2.597592  8.853445  5.685248  4.164969 -2.981424  2.372588 -2.453734  2.737896 
    ##       122       123       125       126       127       129       131       140 
    ##  5.720913 -2.597592  4.344335  2.448979  7.227285  6.872714  5.611772 21.943564 
    ##       171       172       174       175       183       186       211       215 
    ##  9.558049  9.127082  3.194823  2.547771  3.980982  7.879465 -2.279918  7.879465

    atipicos_CUADRA_tecdatos <- as.numeric(names(boxplot(residuos_ROB_CUADRA_tecdatos)$out))

![](README_files/figure-markdown_strict/robust%20est%20FLM-2.png)

    print("The outliers from the robust linear fit are:")

    ## [1] "The outliers from the robust linear fit are:"

    print(atipicos_LINEAL_tecdatos)

    ##  [1]   7  34  35  43  44  45 118 119 121 122 126 127 128 129 130 140 143 168 171
    ## [20] 172 185 186 215

    # 7  34  35  43  44  45 118 119 121 122 126 127 128 129 130 140 143 168 171
    # 172 185 186 215

    print("The common outliers with the functional boxplot are")

    ## [1] "The common outliers with the functional boxplot are"

    intersect(fd1$outpoint[fd1$outpoint <= length(indices_tecdatos)], atipicos_LINEAL_tecdatos)

    ## [1]  35 140

     # 35 140

    print("The outliers from the robust quadratic fit are:")

    ## [1] "The outliers from the robust quadratic fit are:"

    print(atipicos_CUADRA_tecdatos)

    ##  [1]   4   8  10  20  31  34  35  38  40  43  44  45 102 108 117 121 122 123 125
    ## [20] 126 127 129 131 140 171 172 174 175 183 186 211 215

    #  4   8  10  20  31  34  35  38  40  43  44  45 102 108 117 121 122 123 125
    # 126 127 129 131 140 171 172 174 175 183 186 211 215

    print("The common outliers with the functional boxplot are")

    ## [1] "The common outliers with the functional boxplot are"

    intersect(fd1$outpoint[fd1$outpoint <= length(indices_tecdatos)], atipicos_CUADRA_tecdatos)

    ## [1]  35 140

    # 35 140

    ######################################
    # PLOT THE CURVES CORRESPONDING TO 
    # ATYPICAL Residuals from FQM
    #####################################
     
    nombre="atipicos_CUADRA_tecdatos.pdf"
    pdf(nombre,bg='transparent')
    matplot(absorp1$argvals, t(absorp1$data), type = "l", lty = 1,lwd=2,
            xlab = "Wavelength", ylab = "X(t)",col="gray60")

    matplot(absorp1$argvals, t(absorp1$data[atipicos_CUADRA_tecdatos,]), type = "l", 
        lty = 2,lwd=2,
            xlab = "Wavelength", ylab = "X(t)",col="red",add=TRUE)

    dev.off()

    ## png 
    ##   2

    ######################################
    # SAVE THE RESIDUAL BOXPLOT FROM FQM
    #####################################

    nombre="boxplot-rob-quad.pdf"
    pdf(nombre,bg='transparent')

    boxplot(residuos_ROB_CUADRA_tecdatos,col="steelblue")

    dev.off()

    ## png 
    ##   2

    ######################################
    # SAVE THE RESIDUAL BOXPLOT FROM FLM
    #####################################


    nombre="boxplot-rob-lineal.pdf"
    pdf(nombre,bg='transparent')

    boxplot(residuos_ROB_LINEAL_tecdatos,col="steelblue")

    dev.off()

    ## png 
    ##   2

    ##############################################
    ## PLOT BETA FROM FLM FIT
    ############################################## 


    par(mar = c(5,5,2,2))
    plot(tecdatos$t, beta_rob_LINEAL, type = "l", lwd=3,xlab = "Wavelength", pch = 10,
         col = "steelblue", ylab = expression(hat(beta)))

![](README_files/figure-markdown_strict/robust%20est%20FLM-3.png)

    ##############################################
    ## PLOT BETA AND GAMMA FROM FQM FIT
    ############################################## 

    par(mar = c(5,5,2,2))
    plot(tecdatos$t, beta_rob_CUADRA  , type = "l", lwd=3,xlab = "Wavelength", pch = 10,
         col = "steelblue", ylab = expression(hat(beta)))

![](README_files/figure-markdown_strict/robust%20est%20FLM-4.png)

    ene.col=15
    colores=rainbow(ene.col, s = 1, v = 1, start = 0, end = max(1, ene.col - 1)/ene.col, alpha = 1)
     
    eme= min(gamma_rob_CUADRA)
    EME= max(gamma_rob_CUADRA)

    floor(eme)

    ## [1] -403

    donde=seq(eme , EME,length=ene.col)

    wireframe(gamma_rob_CUADRA,scales = list(arrows = FALSE),   
        xlab="t",ylab="s",zlab=" ",   cex=1.2,row.values=tecdatos$t,
        column.values=tecdatos$t,drape = TRUE,aspect = c(1,1) , 
        colorkey = TRUE,at=donde, col.regions=colores)

![](README_files/figure-markdown_strict/robust%20est%20FLM-5.png)

    ##############################################
    ## PLOT residuals vs predicted
    ############################################## 




    par(mar=c(5,4,4,5)+.1)

    bb_L=boxplot(residuos_ROB_LINEAL_tecdatos)

![](README_files/figure-markdown_strict/robust%20est%20FLM-6.png)

    bb_C=boxplot(residuos_ROB_CUADRA_tecdatos)

![](README_files/figure-markdown_strict/robust%20est%20FLM-7.png)

    plot( predichos_ROB_CUADRA, residuos_ROB_CUADRA_tecdatos,
    xlab=expression(hat(y[i])), ylab=expression(r[i]), cex=1.3, cex.lab=1.3)
     abline(h=bb_C$stats[1,1],col="red")
    abline(h=bb_C$stats[5,1],col="red")

![](README_files/figure-markdown_strict/robust%20est%20FLM-8.png)

    par(mar=c(5,4,4,5)+.1)

    plot( predichos_ROB_LINEAL, residuos_ROB_LINEAL_tecdatos, xlab=expression(hat(y[i])), ylab=expression(r[i]), cex=1.3, cex.lab=1.3)
     abline(h=bb_L$stats[1,1],col="red")
    abline(h=bb_L$stats[5,1],col="red")

![](README_files/figure-markdown_strict/robust%20est%20FLM-9.png)

    ########################
    # SAVE THE PLOTS
    ##########################

    nombre="residuo-vs-yhat-rob-lineal.pdf" 
    pdf(nombre,bg='transparent')

    par(mar=c(5,4,4,5)+.1)

    plot( predichos_ROB_LINEAL, residuos_ROB_LINEAL_tecdatos, xlab=expression(hat(y[i])), ylab=expression(r[i]), cex=1.3, cex.lab=1.3)

    dev.off()

    ## png 
    ##   2

    nombre="residuo-vs-yhat-rob-lineal-zoom.pdf" 
    pdf(nombre,bg='transparent')

    par(mar=c(5,4,4,5)+.1)

    plot( predichos_ROB_LINEAL, residuos_ROB_LINEAL_tecdatos, 
    ylim=c(bb_L$stats[1,], bb_L$stats[5,]),xlab=expression(hat(y[i])),
     xlim=c(0,50),ylab=expression(r[i]), cex=1.3, cex.lab=1.3)

    dev.off()

    ## png 
    ##   2

    nombre="residuo-vs-yhat-rob-QUAD.pdf" 
    pdf(nombre,bg='transparent')

    par(mar=c(5,4,4,5)+.1)

    plot( predichos_ROB_CUADRA, residuos_ROB_CUADRA_tecdatos,
    xlab=expression(hat(y[i])), ylab=expression(r[i]), cex=1.3, cex.lab=1.3)

    dev.off()

    ## png 
    ##   2
