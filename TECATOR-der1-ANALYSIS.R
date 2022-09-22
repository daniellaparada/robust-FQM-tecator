###########################################################################
## ---------------------------------------------------------------------
## Datos reales - TECATOR
## http://lib.stat.cmu.edu/datasets/tecator
##
## 20/06/2022
##  ---------------------------------------------------------------------
############################################################################

######################################################################
# ANALISIS CON TODOS LOS DATOS Y DERIVADA PRIMERA
# PARA COMPARAR CLASICO Y ROBUSTO
######################################################################

## Paquetes
library('fda')          
library('robustbase') 		   # lmrob 
library(fda.usc)
library(gdata)                   # uppertriangle
library(lattice)			   # Plots

source('funciones-auxiliares.R') 
source('descomponer3.R')         # covariance decomposition
source('minimizar.R')            # minimization menu
     


## ---------------------------------------------------------------------
## Data and plots
## ---------------------------------------------------------------------

#####################
## Load data
#####################

datos <- data(tecator)
absorp <- tecator$absorp.fdata
absorp1 <- fdata.deriv(absorp,nderiv = 1)  #  computes the first derivative
  
#####################
## Plot the data
#####################

plot(absorp1)  #  first derivative 

matplot(absorp1$argvals, t(absorp1$data), type = "l", lty = 1,
        xlab = "Wavelength", ylab = "X(t)") 

#######################
## Functional Boxplots 
#######################

#######################
# First derivative
#######################

fd1 <- fbplot(t(absorp1$data), x = absorp1$argvals, method =  "MBD",
       xlab = "Wavelength", ylab = "d(Absorbance, 1)",
       xlim = c(850, 1050), ylim = c(-0.02, 0.05))

#######################
# Absorbances
#######################

fd0<-fbplot(t(absorp$data), x = absorp$argvals, method = "MBD",
       xlab = "Wavelength", ylab = "Absorbance",
       xlim = c(850, 1050), ylim = c(2, 5.5))

#######################
# SAVE THE PLOTS
#######################

nombre="tecator-der1.pdf"
pdf(nombre,bg='transparent')

matplot(absorp1$argvals, t(absorp1$data), type = "l", lty = 1,lwd=2,
        xlab = "Wavelength", ylab = "X(t)")
dev.off()



nombre="tecator-der0.pdf"
pdf(nombre,bg='transparent')

matplot(absorp$argvals, t(absorp$data), type = "l", lty = 1,lwd=2,
        xlab = "Wavelength", ylab = "Absorbance")
dev.off()

#######################################################
##  Functional linear model (FLM)
##  fat= alfa_0+<beta, absorp1> + error
#######################################################

set.seed(124) 
covariable <- 'd1'   
 
indices <- 1:215  
 
#######################################
## DATA SET
######################################

tecdatos <- submuestra(indices, covariable = covariable)
 
indices_tecdatos <- indices
dt_tecdatos <- tecdatos$t[2]-tecdatos$t[1]
 


###############################################
## Choose the number of principal directions
## equal to 4, then varexp should be 1
###############################################


freq=4
varexp=1

###############################################
## ROBUST ESTIMATOR UNDER A FLM 
## fat= alfa_0 + < absorp1, beta_0 > + epsilon
###############################################

est_rob_LINEAL <- estimar(y= tecdatos$y,xcenter= tecdatos$x,  
			ttt= tecdatos$t, ajuste='lineal', freq=freq,  cov_type='gerS',  
			fLoss='lmrob',cterho=3.443689, nresamp=5000, varexp=varexp) 

#######################################
# PERCENTAGE EXPLAINED By 4 directions
#######################################

est_rob_LINEAL$porcentaje

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
atipicos_LINEAL_tecdatos <- as.numeric(names(boxplot(residuos_ROB_LINEAL_tecdatos)$out))


names(residuos_ROB_CUADRA_tecdatos) <- 1:length(residuos_ROB_CUADRA_tecdatos)
boxplot(residuos_ROB_CUADRA_tecdatos)$out
atipicos_CUADRA_tecdatos <- as.numeric(names(boxplot(residuos_ROB_CUADRA_tecdatos)$out))

print("The outliers from the robust linear fit are:")
print(atipicos_LINEAL_tecdatos)
# 7  34  35  43  44  45 118 119 121 122 126 127 128 129 130 140 143 168 171
# 172 185 186 215

print("The common outliers with the functional boxplot are")
intersect(fd1$outpoint[fd1$outpoint <= length(indices_tecdatos)], atipicos_LINEAL_tecdatos)
 # 35 140

print("The outliers from the robust quadratic fit are:")
print(atipicos_CUADRA_tecdatos)
#  4   8  10  20  31  34  35  38  40  43  44  45 102 108 117 121 122 123 125
# 126 127 129 131 140 171 172 174 175 183 186 211 215

print("The common outliers with the functional boxplot are")
intersect(fd1$outpoint[fd1$outpoint <= length(indices_tecdatos)], atipicos_CUADRA_tecdatos)
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


######################################
# SAVE THE RESIDUAL BOXPLOT FROM FQM
#####################################

nombre="boxplot-rob-quad.pdf"
pdf(nombre,bg='transparent')

boxplot(residuos_ROB_CUADRA_tecdatos,col="steelblue")

dev.off()


######################################
# SAVE THE RESIDUAL BOXPLOT FROM FLM
#####################################


nombre="boxplot-rob-lineal.pdf"
pdf(nombre,bg='transparent')

boxplot(residuos_ROB_LINEAL_tecdatos,col="steelblue")

dev.off()


##############################################
## PLOT BETA FROM FLM FIT
############################################## 


par(mar = c(5,5,2,2))
plot(tecdatos$t, beta_rob_LINEAL, type = "l", lwd=3,xlab = "Wavelength", pch = 10,
     col = "steelblue", ylab = expression(hat(beta)))

 

##############################################
## PLOT BETA AND GAMMA FROM FQM FIT
############################################## 

par(mar = c(5,5,2,2))
plot(tecdatos$t, beta_rob_CUADRA  , type = "l", lwd=3,xlab = "Wavelength", pch = 10,
     col = "steelblue", ylab = expression(hat(beta)))
 

ene.col=15
colores=rainbow(ene.col, s = 1, v = 1, start = 0, end = max(1, ene.col - 1)/ene.col, alpha = 1)
 
eme= min(gamma_rob_CUADRA)
EME= max(gamma_rob_CUADRA)

floor(eme)

donde=seq(eme , EME,length=ene.col)

wireframe(gamma_rob_CUADRA,scales = list(arrows = FALSE),   
	xlab="t",ylab="s",zlab=" ",   cex=1.2,row.values=tecdatos$t,
	column.values=tecdatos$t,drape = TRUE,aspect = c(1,1) , 
  	colorkey = TRUE,at=donde, col.regions=colores)
  
 
##############################################
## PLOT residuals vs predicted
############################################## 




par(mar=c(5,4,4,5)+.1)

bb_L=boxplot(residuos_ROB_LINEAL_tecdatos)
bb_C=boxplot(residuos_ROB_CUADRA_tecdatos)
plot( predichos_ROB_CUADRA, residuos_ROB_CUADRA_tecdatos,
xlab=expression(hat(y[i])), ylab=expression(r[i]), cex=1.3, cex.lab=1.3)
 abline(h=bb_C$stats[1,1],col="red")
abline(h=bb_C$stats[5,1],col="red")


par(mar=c(5,4,4,5)+.1)

plot( predichos_ROB_LINEAL, residuos_ROB_LINEAL_tecdatos, xlab=expression(hat(y[i])), ylab=expression(r[i]), cex=1.3, cex.lab=1.3)
 abline(h=bb_L$stats[1,1],col="red")
abline(h=bb_L$stats[5,1],col="red")


########################
# SAVE THE PLOTS
##########################

nombre="residuo-vs-yhat-rob-lineal.pdf" 
pdf(nombre,bg='transparent')

par(mar=c(5,4,4,5)+.1)

plot( predichos_ROB_LINEAL, residuos_ROB_LINEAL_tecdatos, xlab=expression(hat(y[i])), ylab=expression(r[i]), cex=1.3, cex.lab=1.3)

dev.off()

nombre="residuo-vs-yhat-rob-lineal-zoom.pdf" 
pdf(nombre,bg='transparent')

par(mar=c(5,4,4,5)+.1)

plot( predichos_ROB_LINEAL, residuos_ROB_LINEAL_tecdatos, 
ylim=c(bb_L$stats[1,], bb_L$stats[5,]),xlab=expression(hat(y[i])),
 xlim=c(0,50),ylab=expression(r[i]), cex=1.3, cex.lab=1.3)

dev.off()



nombre="residuo-vs-yhat-rob-QUAD.pdf" 
pdf(nombre,bg='transparent')

par(mar=c(5,4,4,5)+.1)

plot( predichos_ROB_CUADRA, residuos_ROB_CUADRA_tecdatos,
xlab=expression(hat(y[i])), ylab=expression(r[i]), cex=1.3, cex.lab=1.3)

dev.off()


 


#######################################
## LEAST SQUARE ESTIMATE UNDER FLM
#######################################

est_CL_LINEAL <- estimar(y= tecdatos$y,xcenter= tecdatos$x,  ttt= tecdatos$t, 
			ajuste='lineal', freq=freq,  cov_type='cl',  fLoss='ls',cterho=3.443689, 
			nresamp=5000, varexp=varexp) 



#######################################
# PERCENTAGE EXPLAINED By 4 directions
#######################################

est_CL_LINEAL$porcentaje
#0.9791543


#############################################
# Store the estimates and compute predictions
# for the residuals
#############################################


beta_CL_LINEAL <- est_CL_LINEAL$beta 
phies_CL_LINEAL <- est_CL_LINEAL$autofun

coef_CL_LINEAL <-   tecdatos$x %*% phies_CL_LINEAL * dt_tecdatos

predichos_CL_LINEAL <- est_CL_LINEAL$alfa + coef_CL_LINEAL%*% est_CL_LINEAL$slope_coef 


#######################################
## Estimador  CLASICO  CUADRATICO
#######################################

est_CL_CUADRA <- estimar(y= tecdatos$y,xcenter= tecdatos$x,  ttt= tecdatos$t, 
			ajuste='quadra', freq=freq,  cov_type='cl',  
			fLoss='ls',cterho=3.443689, nresamp=5000, varexp=varexp) 


#############################################
# Store the estimates and compute predictions
# for the residuals
#############################################


beta_CL_CUADRA <- est_CL_CUADRA$beta
gamma_CL_CUADRA <- est_CL_CUADRA$gamma
phies_CL <- est_CL_CUADRA$autofun
 
coef_CL_CUADRA <-   tecdatos$x %*% phies_CL * dt_tecdatos

predichos_CL_CUADRA <- est_CL_CUADRA$alfa + coef_CL_CUADRA%*% est_CL_CUADRA$slope_coef + diag(coef_CL_CUADRA%*% est_CL_CUADRA$gama_coef%*% t(coef_CL_CUADRA))



#############################################
#  RESIDUALS FROM THE LEAST SQUARES FIT
#############################################
 
 
residuos_CL_CUADRA_tecdatos <- tecdatos$y - predichos_CL_CUADRA

residuos_CL_LINEAL_tecdatos <- tecdatos$y - predichos_CL_LINEAL



##############################################
## PLOTS OF BETA UNDER FLM
############################################## 

par(mar = c(5,5,2,2))
plot(tecdatos$t, beta_CL_LINEAL, type = "l", lwd=3,xlab = "Wavelength", pch = 10,
     col = "red", ylab = expression(hat(beta)),
	 ylim=c(min(c(beta_CL_LINEAL,beta_rob_LINEAL)), max(c(beta_CL_LINEAL,beta_rob_LINEAL))))
 


##############################################
## COMPARE WITH THE ROBUST ONE
############################################## 

par(mar = c(5,5,2,2))
plot(tecdatos$t, beta_CL_LINEAL, type = "l", lwd=3,xlab = "Wavelength", pch = 10,
     col = "red", ylab = expression(hat(beta)),
	 ylim=c(min(c(beta_CL_LINEAL,beta_rob_LINEAL)), max(c(beta_CL_LINEAL,beta_rob_LINEAL))))
 

lines(tecdatos$t, beta_rob_LINEAL, type = "l", lwd=3,xlab = "Wavelength", pch = 10,
     col = "steelblue")



##############################################
## PLOTS OF BETA AND GAMMA UNDER FQM
############################################## 

par(mar = c(5,5,2,2))
plot(tecdatos$t, beta_CL_CUADRA, type = "l", lwd=3,xlab = "Wavelength", pch = 10,
     col = "tomato", ylab = expression(hat(beta)),
     ylim=c(min(c(beta_CL_CUADRA,beta_rob_CUADRA)), max(c(beta_CL_CUADRA,beta_rob_CUADRA))))
 

####################################################
## COMPARE ESTIMATES OF BETA LS WITH THE ROBUST ONE
####################################################

par(mar = c(5,5,2,2))
plot(tecdatos$t, beta_CL_CUADRA, type = "l", lwd=3,xlab = "Wavelength", pch = 10,
     col = "tomato", ylab = expression(hat(beta)),
     ylim=c(min(c(beta_CL_CUADRA,beta_rob_CUADRA)), max(c(beta_CL_CUADRA,beta_rob_CUADRA))))
 

lines(tecdatos$t, beta_rob_CUADRA, type = "l", lwd=3,xlab = "Wavelength", pch = 10,
     col = "steelblue")
 

#####################################
# PLOT OF GAMMA
####################################

ene.col=15
colores=rainbow(ene.col, s = 1, v = 1, start = 0, end = max(1, ene.col - 1)/ene.col, alpha = 1)
 
eme= min(gamma_CL_CUADRA)
eMe= max(gamma_CL_CUADRA)

floor(eme)

donde=seq(eme ,eMe  ,length=ene.col)
 
wireframe(gamma_CL_CUADRA,scales = list(arrows = FALSE),   
		xlab="t",ylab="s",zlab=" ",   cex=1.2,
		row.values=tecdatos$t,column.values=tecdatos$t,drape = TRUE,
  		colorkey = TRUE,at=donde, col.regions=colores)
   


#####################################
# PLOT Residuals versus predicted  
# Using classical procedure
####################################

par(mar=c(5,4,4,5)+.1)

plot( predichos_CL_CUADRA, residuos_CL_CUADRA_tecdatos,
xlab=expression(hat(y[i])), ylab=expression(r[i]), cex=1.3, cex.lab=1.3)


par(mar=c(5,4,4,5)+.1)

plot( predichos_CL_LINEAL,residuos_CL_LINEAL_tecdatos, xlab=expression(hat(y[i])), ylab=expression(r[i]), cex=1.3, cex.lab=1.3)



########################
# SAVE PLOTS
##########################

nombre="residuo-vs-yhat-CL-lineal.pdf" 
pdf(nombre,bg='transparent')

par(mar=c(5,4,4,5)+.1)

plot( predichos_CL_LINEAL,residuos_CL_LINEAL_tecdatos, xlab=expression(hat(y[i])), ylab=expression(r[i]), cex=1.3, cex.lab=1.3)


dev.off()


nombre="residuo-vs-yhat-CL-QUAD.pdf" 
pdf(nombre,bg='transparent')

par(mar=c(5,4,4,5)+.1)

plot( predichos_CL_CUADRA,residuos_CL_CUADRA_tecdatos, xlab=expression(hat(y[i])), ylab=expression(r[i]), cex=1.3, cex.lab=1.3)


dev.off()

###################################
# SAVE PLOTS of estimates of beta
###################################

####FQM

nombre="beta-CL-ROB-QUAD.pdf" 
pdf(nombre,bg='transparent')

par(mar = c(5,5,2,2))
plot(tecdatos$t, beta_CL_CUADRA, type = "l", lwd=3,xlab = "Wavelength", pch = 10,
     col = "tomato", ylab = expression(hat(beta)),
     ylim=c(min(c(beta_CL_CUADRA,beta_rob_CUADRA)), max(c(beta_CL_CUADRA,beta_rob_CUADRA))))
 
lines(tecdatos$t, beta_rob_CUADRA, type = "l", lwd=3,xlab = "Wavelength", pch = 10,
     col = "blue")
dev.off()

####FLM

nombre="beta-CL-ROB-LINEAL.pdf" 
pdf(nombre,bg='transparent')

par(mar = c(5,5,2,2))
plot(tecdatos$t, beta_CL_LINEAL, type = "l", lwd=3,xlab = "Wavelength", pch = 10,
     col = "red", ylab = expression(hat(beta)),
	 ylim=c(min(c(beta_CL_LINEAL,beta_rob_LINEAL)), max(c(beta_CL_LINEAL,beta_rob_LINEAL))))
 
lines(tecdatos$t, beta_rob_LINEAL, type = "l", lwd=3,xlab = "Wavelength", pch = 10,
     col = "blue")

dev.off()

############################################################# 
# REPEAT THE ANALYSIS FOR THE CLASSICAL PROCEDURE UNDER FQM
# WITHOUT THE OUTLIERS DETECTED  BY THE ROBUST METHOD
#############################################################

############################ 
## Data without outliers
############################

tecdatos_sin_out <- submuestra(indices_tecdatos[-atipicos_CUADRA_tecdatos], covariable = covariable)
 


#######################################
## CLASSICAL ESTIMATOR WITHOUT OULTIERS
## FLM
#######################################

est_CL_LINEAL_sin_out <- estimar(y= tecdatos_sin_out$y,xcenter= tecdatos_sin_out$x,  ttt= tecdatos_sin_out$t, 
				ajuste='lineal', freq=freq,  cov_type='cl',  fLoss='ls',cterho=3.443689, 
				nresamp=5000, varexp=varexp) 

est_CL_LINEAL_sin_out$porcentaje
#0.9774696

#####################################################
# Store parameters and compute predicted values
######################################################


beta_CL_LINEAL_sin_out <- est_CL_LINEAL_sin_out$beta 
phies_CL_LINEAL_sin_out <- est_CL_LINEAL_sin_out$autofun

coef_CL_LINEAL_sin_out <-   tecdatos_sin_out$x %*% phies_CL_LINEAL_sin_out * dt_tecdatos

predichos_CL_LINEAL_sin_out <- est_CL_LINEAL_sin_out$alfa + coef_CL_LINEAL_sin_out%*% est_CL_LINEAL_sin_out$slope_coef 


######################################
## CLASSICAL ESTIMATOR WITHOUT OULTIERS
## FQM
#######################################

est_CL_CUADRA_sin_out <- estimar(y= tecdatos_sin_out$y,xcenter= tecdatos_sin_out$x,  
				ttt= tecdatos_sin_out$t, 
				ajuste='quadra', freq=freq,  cov_type='cl',  
				fLoss='ls',cterho=3.443689, nresamp=5000, varexp=varexp) 

######################################################
# Store parameters and compute predicted values
######################################################

beta_CL_CUADRA_sin_out <- est_CL_CUADRA_sin_out$beta
gamma_CL_CUADRA_sin_out <- est_CL_CUADRA_sin_out$gamma
phies_CL_sin_out <- est_CL_CUADRA_sin_out$autofun
 
coef_CL_CUADRA_sin_out <-   tecdatos_sin_out$x %*% phies_CL_sin_out * dt_tecdatos

predichos_CL_CUADRA_sin_out <- est_CL_CUADRA_sin_out$alfa + coef_CL_CUADRA_sin_out%*% est_CL_CUADRA_sin_out$slope_coef + 
diag(coef_CL_CUADRA_sin_out%*% est_CL_CUADRA_sin_out$gama_coef%*% t(coef_CL_CUADRA_sin_out))


######################################################
## RESIDUALS FROM THE CLASSICAL FIT WITHOUT OUTLIERS
###################################################### 

residuos_CL_CUADRA_tecdatos_sin_out <- tecdatos_sin_out$y - predichos_CL_CUADRA_sin_out

residuos_CL_LINEAL_tecdatos_sin_out <- tecdatos_sin_out$y - predichos_CL_LINEAL_sin_out


####################################
# PLOTS
######################################
par(mar=c(5,4,4,5)+.1)

plot( predichos_CL_CUADRA_sin_out, residuos_CL_CUADRA_tecdatos_sin_out,
xlab=expression(hat(y[i])), ylab=expression(r[i]), cex=1.3, cex.lab=1.3)


par(mar=c(5,4,4,5)+.1)

plot( predichos_CL_LINEAL_sin_out,residuos_CL_LINEAL_tecdatos_sin_out, xlab=expression(hat(y[i])), ylab=expression(r[i]), cex=1.3, cex.lab=1.3)



par(mar = c(5,5,2,2))
plot(tecdatos$t, beta_CL_CUADRA, type = "l",  lwd=3,xlab = "Wavelength", pch = 10,
     col = "tomato", ylab = expression(hat(beta)),
     ylim=c(min(c(beta_CL_CUADRA,beta_rob_CUADRA)), max(c(beta_CL_CUADRA,beta_rob_CUADRA))))
 
lines(tecdatos$t, beta_rob_CUADRA, type = "l", lwd=3,xlab = "Wavelength", pch = 10,
     col = "blue")
 
lines(tecdatos$t, beta_CL_CUADRA_sin_out, type = "l", lty=2,lwd=4,xlab = "Wavelength", pch = 10,
     col = "pink")

par(mar = c(5,5,2,2))
plot(tecdatos$t, beta_CL_LINEAL, type = "l", lwd=3,xlab = "Wavelength", pch = 10,
     col = "red", ylab = expression(hat(beta)),
	 ylim=c(min(c(beta_CL_LINEAL,beta_rob_LINEAL)), max(c(beta_CL_LINEAL,beta_rob_LINEAL))))
 
lines(tecdatos$t, beta_rob_LINEAL, type = "l", lwd=3,xlab = "Wavelength", pch = 10,
     col = "blue")

lines(tecdatos$t, beta_CL_LINEAL_sin_out, type = "l", lty=2,lwd=4,xlab = "Wavelength", pch = 10,
     col = "pink")






ene.col=15
colores=rainbow(ene.col, s = 1, v = 1, start = 0, end = max(1, ene.col - 1)/ene.col, alpha = 1)
 
eme= min(c(gamma_CL_CUADRA_sin_out,gamma_CL_CUADRA,gamma_rob_CUADRA) )
eMe= max(c(gamma_CL_CUADRA_sin_out,gamma_CL_CUADRA,gamma_rob_CUADRA) )

floor(eme)

donde=seq(eme ,eMe  ,length=ene.col)

wireframe(gamma_CL_CUADRA,scales = list(arrows = FALSE),   
		xlab="t",ylab="s",zlab=" ",   cex=1.2, zlim=c(eme ,eMe),screen=list(z= 30,x=- 70),
		row.values=tecdatos$t,column.values=tecdatos$t,drape = TRUE,
  		colorkey = TRUE,at=donde)

windows()

wireframe(gamma_rob_CUADRA ,scales = list(arrows = FALSE),   
		xlab="t",ylab="s",zlab=" ",   cex=1.2,zlim=c(eme ,eMe),screen=list(z= 30,x=- 70),
		row.values=tecdatos$t,column.values=tecdatos$t,drape = TRUE,
  		colorkey = TRUE, at=donde)

windows()


wireframe(gamma_CL_CUADRA_sin_out,scales = list(arrows = FALSE),   
		xlab="t",ylab="s",zlab=" ",   cex=1.2,zlim=c(eme ,eMe),screen=list(z= 30,x=- 70),
		row.values=tecdatos$t,column.values=tecdatos$t,drape = TRUE,
  		colorkey = TRUE, at=donde)


 


#######################################
# SURFACE DIFERENCES
##########################################

 

###################################
#SURFACES
###################################
 

DIF_CL_ROB <- gamma_CL_CUADRA - gamma_rob_CUADRA
DIF_CL_sin_out_ROB <- gamma_CL_CUADRA_sin_out- gamma_rob_CUADRA

eme_D=min(c(DIF_CL_ROB,DIF_CL_sin_out_ROB))
EME_D=max(c(DIF_CL_ROB,DIF_CL_sin_out_ROB))



donde_D=seq(eme_D ,EME_D  ,length=ene.col)

wireframe(DIF_CL_ROB ,scales = list(arrows = FALSE),   
	xlab="t",ylab="s",zlab=" ",   cex=1.2,zlim=c(eme_D ,EME_D),
	row.values=tecdatos$t,column.values=tecdatos$t,drape = TRUE,
  	colorkey = TRUE,aspect = c(1,1),screen=list(z= 30,x=- 70),
 	at=donde_D, col.regions=colores) 

  

windows()

wireframe(DIF_CL_sin_out_ROB ,scales = list(arrows = FALSE),   
	xlab="t",ylab="s",zlab=" ",   cex=1.2,zlim=c(eme_D ,EME_D),
	row.values=tecdatos$t,column.values=tecdatos$t,drape = TRUE,
  	colorkey = TRUE, aspect = c(1,1),screen=list(z= 30,x=- 70),
	at=donde_D, col.regions=colores) 
 

 


########################
# SAVE THE PLOTS AS PDF
##########################

nombre="residuo-vs-yhat-CL-sin-out-lineal.pdf" 
pdf(nombre,bg='transparent')
 

par(mar=c(5,4,4,5)+.1)

plot( predichos_CL_LINEAL_sin_out,residuos_CL_LINEAL_tecdatos_sin_out, xlab=expression(hat(y[i])), ylab=expression(r[i]), cex=1.3, cex.lab=1.3)

dev.off()

nombre="residuo-vs-yhat-CL-sin-out-QUAD.pdf" 
pdf(nombre,bg='transparent')
 
par(mar=c(5,4,4,5)+.1)

plot( predichos_CL_CUADRA_sin_out, residuos_CL_CUADRA_tecdatos_sin_out,
xlab=expression(hat(y[i])), ylab=expression(r[i]), cex=1.3, cex.lab=1.3,
ylim=c(-2,2))

dev.off()


nombre="residuo-vs-yhat-rob-QUAD-zoom.pdf" 
pdf(nombre,bg='transparent')
 
par(mar=c(5,4,4,5)+.1)

plot( predichos_CL_CUADRA_sin_out, residuos_CL_CUADRA_tecdatos_sin_out,
xlab=expression(hat(y[i])), ylab=expression(r[i]), cex=1.3, cex.lab=1.3,
ylim=c(-2,2))
#abline(h=bb_C$stats[1,1],col="red")
#abline(h=bb_C$stats[5,1],col="red")

dev.off()


nombre="beta-CL-CL-SO-ROB-QUAD.pdf"
pdf(nombre,bg='transparent')

par(mar = c(5,5,2,2))
plot(tecdatos$t, beta_CL_CUADRA, type = "l",  lwd=3,xlab = "Wavelength", pch = 10,
     col = "tomato", ylab = expression(hat(beta)),
     ylim=c(min(c(beta_CL_CUADRA,beta_rob_CUADRA)), max(c(beta_CL_CUADRA,beta_rob_CUADRA))))
 
lines(tecdatos$t, beta_rob_CUADRA, type = "l", lwd=3,xlab = "Wavelength", pch = 10,
     col = "blue")
 
lines(tecdatos$t, beta_CL_CUADRA_sin_out, type = "l", lty=2,lwd=4,xlab = "Wavelength", pch = 10,
     col = "pink2")

dev.off()


   


nombre="beta-CL-ROB-QUAD-new2.pdf"
pdf(nombre,bg='transparent')

par(mar = c(5,5,2,2))
plot(tecdatos$t, beta_CL_CUADRA, type = "l",  lwd=3,xlab = "Wavelength", pch = 10,
     col = "red", ylab = expression(hat(beta)),
     ylim=c(min(c(beta_CL_CUADRA,beta_rob_CUADRA)), max(c(beta_CL_CUADRA,beta_rob_CUADRA))))
 
lines(tecdatos$t, beta_rob_CUADRA, type = "l", lwd=3,xlab = "Wavelength", pch = 10,
     col = "steelblue")
  

dev.off()

nombre="beta-CL-CL-SO-ROB-LINEAL.pdf"
pdf(nombre,bg='transparent')



par(mar = c(5,5,2,2))
plot(tecdatos$t, beta_CL_LINEAL, type = "l", lwd=3,xlab = "Wavelength", pch = 10,
     col = "red", ylab = expression(hat(beta)),
	 ylim=c(min(c(beta_CL_LINEAL,beta_rob_LINEAL)), max(c(beta_CL_LINEAL,beta_rob_LINEAL))))
 
lines(tecdatos$t, beta_rob_LINEAL, type = "l", lwd=3,xlab = "Wavelength", pch = 10,
     col = "blue")

lines(tecdatos$t, beta_CL_LINEAL_sin_out, type = "l", lty=2,lwd=4,xlab = "Wavelength", pch = 10,
     col = "pink2")


dev.off()


 

nombre="beta-CL-ROB-LINEAL-new2.pdf"
pdf(nombre,bg='transparent')



par(mar = c(5,5,2,2))
plot(tecdatos$t, beta_CL_LINEAL, type = "l", lwd=3,xlab = "Wavelength", pch = 10,
     col = "red", ylab = expression(hat(beta)),
	 ylim=c(min(c(beta_CL_LINEAL,beta_rob_LINEAL)), max(c(beta_CL_LINEAL,beta_rob_LINEAL))))
 
lines(tecdatos$t, beta_rob_LINEAL, type = "l", lwd=3,xlab = "Wavelength", pch = 10,
     col = "steelblue")
 


dev.off()


###################################
#DIFERENCIAS COLORES
###################################
 
 

nombre="DIF-GAMMA-CL-ROB.pdf"
pdf(nombre,bg='transparent')

wireframe(DIF_CL_ROB ,scales = list(arrows = FALSE),   
	xlab="t",ylab="s",zlab=" ",   cex=1.2,zlim=c(eme_D ,EME_D),
	row.values=tecdatos$t,column.values=tecdatos$t,drape = TRUE,
  	colorkey = TRUE,aspect = c(1,1),screen=list(z= 30,x=- 70),
 	at=donde_D, col.regions=colores) 

 
dev.off()
 


nombre="DIF-GAMMA-CL-SO-ROB.pdf"
pdf(nombre,bg='transparent')


wireframe(DIF_CL_sin_out_ROB ,scales = list(arrows = FALSE),   
	xlab="t",ylab="s",zlab=" ",   cex=1.2,zlim=c(eme_D ,EME_D),
	row.values=tecdatos$t,column.values=tecdatos$t,drape = TRUE,
  	colorkey = TRUE, aspect = c(1,1),screen=list(z= 30,x=- 70),
	at=donde_D, col.regions=colores) 


dev.off()

eme1=min(DIF_CL_sin_out_ROB)
EME1= max(DIF_CL_sin_out_ROB)


donde1=seq(eme1 ,EME1  ,length=ene.col)



nombre="DIF-GAMMA-CL-SO-ROB-zoom.pdf"
pdf(nombre,bg='transparent')

wireframe(DIF_CL_sin_out_ROB ,scales = list(arrows = FALSE),   
	xlab="t",ylab="s",zlab=" ",   cex=1.2,zlim=c(eme1 ,EME1),
	row.values=tecdatos$t,column.values=tecdatos$t,drape = TRUE,
  	colorkey = TRUE, aspect = c(1,1),screen=list(z= 30,x=- 70),
	at=donde1, col.regions=colores) 

dev.off()


#####################################
# GAMMA EN COLORES
######################################



nombre="GAMMA-CL.pdf"
pdf(nombre,bg='transparent')


wireframe(gamma_CL_CUADRA,scales = list(arrows = FALSE),   
	xlab="t",ylab="s",zlab=" ",   cex=1.2,zlim=c(eme ,eMe),
	row.values=tecdatos$t,column.values=tecdatos$t,drape = TRUE,
  	colorkey = TRUE,aspect = c(1,1),screen=list(z= 30,x=- 70),
 	at=donde, col.regions=colores) 
dev.off()
  
 

nombre="GAMMA-ROB.pdf"
pdf(nombre,bg='transparent')

wireframe(gamma_rob_CUADRA ,scales = list(arrows = FALSE),   
	xlab="t",ylab="s",zlab=" ",   cex=1.2,zlim=c(eme ,eMe),
	row.values=tecdatos$t,column.values=tecdatos$t,drape = TRUE,
  	colorkey = TRUE, aspect = c(1,1),screen=list(z= 30,x=- 70),
	at=donde, col.regions=colores) 

dev.off()


nombre="GAMMA-CL-SO.pdf"
pdf(nombre,bg='transparent')

wireframe(gamma_CL_CUADRA_sin_out,scales = list(arrows = FALSE),   
	xlab="t",ylab="s",zlab=" ",   cex=1.2,zlim=c(eme ,eMe),
	row.values=tecdatos$t,column.values=tecdatos$t,drape = TRUE,
  	colorkey = TRUE, aspect = c(1,1),screen=list(z= 30,x=- 70),
 	at=donde, col.regions=colores) 

dev.off()

