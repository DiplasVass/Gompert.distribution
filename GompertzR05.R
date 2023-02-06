.libPaths("/home/testuser/Rlib")
install.packages('VGAM')
library(VGAM)

#------------------------------------------------------------------------------------------------
                                    ##Simulation Study##
#------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------
                                       # phi = 0.5 #
#------------------------------------------------------------------------------------------------
phi = 0.5

L =5000 # We want the process to be repeated 5000 times 

t = seq(0.1,0.5,0.1)

Rbty = function(x) { #Construct the Real Reliability Value for known phi

YY = exp(phi*(1-exp(x)))

return(YY) 
			 }

#==========================================================================================
                                        #n=15#
#==========================================================================================
n=15

## Monte Carlo Simulations (L = 5000) for n=15

Dat15 = replicate(L,rgompertz(15, scale = 1, shape=0.5))

#-------------------------------------------------------------------------------------------
                        ##Statistical Estimation for phi##                                        
#-------------------------------------------------------------------------------------------

##### Estimators #####

#-------------------------------------------------------------------------------------------

# We foumd the T estimator which help us for the estimation of phi

T15 = 0

for(i in 1:L) {  # Find the 5000 values of T estimator  

T15[i] = sum(1-exp(Dat15[,i]))

		   }

T15

#-------------------------------------------
                  #MLE#
#-------------------------------------------

Rb = 0

RMLE15  = function(x) {

Rb = exp((-n/T15)*(1-exp(x)))

return(Rb)			   
 			     }


#-----------------
    #IMSE#
#-----------------

klasma1 = 1/L

klasma2 = 1/length(t)

XX = rep(NA,5)

for(j in 1:length(t)) {

XX[j] = sum((RMLE15(t[j])-Rbty(t[j]))^2)/L

                       }

errmle15 = klasma2*sum(XX)
errmle15

#-------------------------------------------
                  #UMVUE#
#-------------------------------------------

Rb = 0

RUMV15  = function(x) {

Rb = exp((-(n-1)/T15)*(1-exp(x)))

return(Rb)			   
 			     }


#-----------------
    #IMSE#
#-----------------

klasma1 = 1/L

klasma2 = 1/length(t)

XX = rep(NA,5)

for(j in 1:length(t)) {

XX[j] = sum((RUMV15(t[j])-Rbty(t[j]))^2)/L

                       }

errumv15 = klasma2*sum(XX)
errumv15

#-------------------------------------------
                  #MSE#
#-------------------------------------------

Rb = 0

RMSE15  = function(x) {

Rb = exp((-(n-2)/T15)*(1-exp(x)))
              
return(Rb)			   
 			     }


#-----------------
    #IMSE#
#-----------------

klasma1 = 1/L

klasma2 = 1/length(t)

XX = rep(NA,5)

for(j in 1:length(t)) {

XX[j] = sum((RMSE15(t[j])-Rbty(t[j]))^2)/L

                       }

errmse15 = klasma2*sum(XX)
errmse15

#-------------------------------------------
                  #Jeffreys#
#-------------------------------------------

Rb = 0

RJEF15  = function(x) {

for(i in 1:L) {

Rb = ((-T15+2*(1-exp(x)))/(-T15+(1-exp(x))))^n

              }
return(Rb)			   
 			     }


#-----------------
    #IMSE#
#-----------------

klasma1 = 1/L

klasma2 = 1/length(t)

XX = rep(NA,5)

for(j in 1:length(t)) {

XX[j] = sum((RJEF15(t[j])-Rbty(t[j]))^2)/L

                       }

errjef15 = klasma2*sum(XX)
errjef15

#-------------------------------------------
                  #Gamma#
#-------------------------------------------

#-------------------------------------------
                 #gama=0.8#
#-------------------------------------------
gama=0.8
#------------------
#delta=0.5#
#------------------

delta=0.5

Rb = 0

RG1D115  = function(x) {

for(i in 1:L) {

Rb = ((delta-T15+2*(1-exp(x)))/(delta-T15+(1-exp(x))))^(n+gama)

              }
return(Rb)			   
 			     }


#-----------------
    #IMSE#
#-----------------

klasma1 = 1/L

klasma2 = 1/length(t)

XX = rep(NA,5)

for(j in 1:length(t)) {

XX[j] = sum((RG1D115(t[j])-Rbty(t[j]))^2)/L

                       }

errg1d115 = klasma2*sum(XX)
errg1d115

#------------------
#delta=3.0#
#------------------

delta=3.0

Rb = 0

RG1D215  = function(x) {

for(i in 1:L) {

Rb = ((delta-T15+2*(1-exp(x)))/(delta-T15+(1-exp(x))))^(n+gama)

              }
return(Rb)			   
 			     }


#-----------------
    #IMSE#
#-----------------

klasma1 = 1/L

klasma2 = 1/length(t)

XX = rep(NA,5)

for(j in 1:length(t)) {

XX[j] = sum((RG1D215(t[j])-Rbty(t[j]))^2)/L

                       }

errg1d215 = klasma2*sum(XX)
errg1d215

#-------------------------------------------
                 #gama=3.0#
#--------------------------------------------
gama=3.0
#------------------
#delta=0.5#
#------------------

delta=0.5

Rb = 0

RG2D115  = function(x) {

for(i in 1:L) {

Rb = ((delta-T15+2*(1-exp(x)))/(delta-T15+(1-exp(x))))^(n+gama)

              }
return(Rb)			   
 			     }


#-----------------
    #IMSE#
#-----------------

klasma1 = 1/L

klasma2 = 1/length(t)

XX = rep(NA,5)

for(j in 1:length(t)) {

XX[j] = sum((RG2D115(t[j])-Rbty(t[j]))^2)/L

                       }

errg2d115 = klasma2*sum(XX)
errg2d115

#------------------
#delta=3.0#
#------------------

delta=3.0

Rb = 0

RG2D215  = function(x) {

for(i in 1:L) {

Rb = ((delta-T15+2*(1-exp(x)))/(delta-T15+(1-exp(x))))^(n+gama)

              }
return(Rb)			   
 			     }


#-----------------
    #IMSE#
#-----------------

klasma1 = 1/L

klasma2 = 1/length(t)

XX = rep(NA,5)

for(j in 1:length(t)) {

XX[j] = sum((RG2D215(t[j])-Rbty(t[j]))^2)/L

                       }

errg2d215 = klasma2*sum(XX)
errg2d215

tb15 = c(errmle15,errumv15,errmse15,errjef15,errg1d115,errg1d215,errg2d115,errg2d215)
tb15

#==========================================================================================
                                        #n=50#
#==========================================================================================
n=50

## Monte Carlo Simulations (L = 5000) for n=50

Dat50 = replicate(L,rgompertz(50, scale = 1, shape=0.5))

#-------------------------------------------------------------------------------------------
                        ##Statistical Estimation for phi##                                        
#-------------------------------------------------------------------------------------------

##### Estimators #####

#-------------------------------------------------------------------------------------------

# We foumd the T estimator which help us for the estimation of phi

T50 = 0

for(i in 1:L) {  # Find the 5000 values of T estimator  

T50[i] = sum(1-exp(Dat50[,i]))

		   }

T50

#-------------------------------------------
                  #MLE#
#-------------------------------------------

Rb = 0

RMLE50  = function(x) {

Rb = exp((-n/T50)*(1-exp(x)))

return(Rb)			   
 			     }


#-----------------
    #IMSE#
#-----------------

klasma1 = 1/L

klasma2 = 1/length(t)

XX = rep(NA,5)

for(j in 1:length(t)) {

XX[j] = sum((RMLE50(t[j])-Rbty(t[j]))^2)/L

                       }

errmle50 = klasma2*sum(XX)
errmle50

#-------------------------------------------
                  #UMVUE#
#-------------------------------------------

Rb = 0

RUMV50  = function(x) {

Rb = exp((-(n-1)/T50)*(1-exp(x)))

return(Rb)			   
 			     }


#-----------------
    #IMSE#
#-----------------

klasma1 = 1/L

klasma2 = 1/length(t)

XX = rep(NA,5)

for(j in 1:length(t)) {

XX[j] = sum((RUMV50(t[j])-Rbty(t[j]))^2)/L

                       }

errumv50 = klasma2*sum(XX)
errumv50

#-------------------------------------------
                  #MSE#
#-------------------------------------------

Rb = 0

RMSE50  = function(x) {

Rb = exp((-(n-2)/T50)*(1-exp(x)))
              
return(Rb)			   
 			     }


#-----------------
    #IMSE#
#-----------------

klasma1 = 1/L

klasma2 = 1/length(t)

XX = rep(NA,5)

for(j in 1:length(t)) {

XX[j] = sum((RMSE50(t[j])-Rbty(t[j]))^2)/L

                       }

errmse50 = klasma2*sum(XX)
errmse50

#-------------------------------------------
                  #Jeffreys#
#-------------------------------------------

Rb = 0

RJEF50  = function(x) {

for(i in 1:L) {

Rb = ((-T50+2*(1-exp(x)))/(-T50+(1-exp(x))))^n

              }
return(Rb)			   
 			     }


#-----------------
    #IMSE#
#-----------------

klasma1 = 1/L

klasma2 = 1/length(t)

XX = rep(NA,5)

for(j in 1:length(t)) {

XX[j] = sum((RJEF50(t[j])-Rbty(t[j]))^2)/L

                       }

errjef50 = klasma2*sum(XX)
errjef50

#-------------------------------------------
                  #Gamma#
#-------------------------------------------

#-------------------------------------------
                 #gama=0.8#
#-------------------------------------------
gama=0.8
#------------------
#delta=0.5#
#------------------

delta=0.5

Rb = 0

RG1D150  = function(x) {

for(i in 1:L) {

Rb = ((delta-T50+2*(1-exp(x)))/(delta-T50+(1-exp(x))))^(n+gama)

              }
return(Rb)			   
 			     }


#-----------------
    #IMSE#
#-----------------

klasma1 = 1/L

klasma2 = 1/length(t)

XX = rep(NA,5)

for(j in 1:length(t)) {

XX[j] = sum((RG1D150(t[j])-Rbty(t[j]))^2)/L

                       }

errg1d150 = klasma2*sum(XX)
errg1d150

#------------------
#delta=3.0#
#------------------

delta=3.0

Rb = 0

RG1D250  = function(x) {

for(i in 1:L) {

Rb = ((delta-T50+2*(1-exp(x)))/(delta-T50+(1-exp(x))))^(n+gama)

              }
return(Rb)			   
 			     }


#-----------------
    #IMSE#
#-----------------

klasma1 = 1/L

klasma2 = 1/length(t)

XX = rep(NA,5)

for(j in 1:length(t)) {

XX[j] = sum((RG1D250(t[j])-Rbty(t[j]))^2)/L

                       }

errg1d250 = klasma2*sum(XX)
errg1d250

#-------------------------------------------
                 #gama=3.0#
#--------------------------------------------
gama=3.0
#------------------
#delta=0.5#
#------------------

delta=0.5

Rb = 0

RG2D150  = function(x) {

for(i in 1:L) {

Rb = ((delta-T50+2*(1-exp(x)))/(delta-T50+(1-exp(x))))^(n+gama)

              }
return(Rb)			   
 			     }


#-----------------
    #IMSE#
#-----------------

klasma1 = 1/L

klasma2 = 1/length(t)

XX = rep(NA,5)

for(j in 1:length(t)) {

XX[j] = sum((RG2D150(t[j])-Rbty(t[j]))^2)/L

                       }

errg2d150 = klasma2*sum(XX)
errg2d150

#------------------
#delta=3.0#
#------------------

delta=3.0

Rb = 0

RG2D250  = function(x) {

for(i in 1:L) {

Rb = ((delta-T50+2*(1-exp(x)))/(delta-T50+(1-exp(x))))^(n+gama)

              }
return(Rb)			   
 			     }


#-----------------
    #IMSE#
#-----------------

klasma1 = 1/L

klasma2 = 1/length(t)

XX = rep(NA,5)

for(j in 1:length(t)) {

XX[j] = sum((RG2D250(t[j])-Rbty(t[j]))^2)/L

                       }

errg2d250 = klasma2*sum(XX)
errg2d250

tb50 = c(errmle50,errumv50,errmse50,errjef50,errg1d150,errg1d250,errg2d150,errg2d250)
tb50

#==========================================================================================
                                        #n=100#
#==========================================================================================
n=100

## Monte Carlo Simulations (L = 5000) for n=50

Dat100 = replicate(L,rgompertz(100, scale = 1, shape=0.5))

#-------------------------------------------------------------------------------------------
                        ##Statistical Estimation for phi##                                        
#-------------------------------------------------------------------------------------------

##### Estimators #####

#-------------------------------------------------------------------------------------------

# We foumd the T estimator which help us for the estimation of phi

T100 = 0

for(i in 1:L) {  # Find the 5000 values of T estimator  

T100[i] = sum(1-exp(Dat100[,i]))

		   }

T100

#-------------------------------------------
                  #MLE#
#-------------------------------------------

Rb = 0

RMLE100  = function(x) {

Rb = exp((-n/T100)*(1-exp(x)))

return(Rb)			   
 			     }


#-----------------
    #IMSE#
#-----------------

klasma1 = 1/L

klasma2 = 1/length(t)

XX = rep(NA,5)

for(j in 1:length(t)) {

XX[j] = sum((RMLE100(t[j])-Rbty(t[j]))^2)/L

                       }

errmle100 = klasma2*sum(XX)
errmle100

#-------------------------------------------
                  #UMVUE#
#-------------------------------------------

Rb = 0

RUMV100  = function(x) {

Rb = exp((-(n-1)/T100)*(1-exp(x)))

return(Rb)			   
 			     }

#-----------------
    #IMSE#
#-----------------

klasma1 = 1/L

klasma2 = 1/length(t)

XX = rep(NA,5)

for(j in 1:length(t)) {

XX[j] = sum((RUMV100(t[j])-Rbty(t[j]))^2)/L

                       }

errumv100 = klasma2*sum(XX)
errumv100

#-------------------------------------------
                  #MSE#
#-------------------------------------------

Rb = 0

RMSE100  = function(x) {

Rb = exp((-(n-2)/T100)*(1-exp(x)))
              
return(Rb)			   
 			     }


#-----------------
    #IMSE#
#-----------------

klasma1 = 1/L

klasma2 = 1/length(t)

XX = rep(NA,5)

for(j in 1:length(t)) {

XX[j] = sum((RMSE100(t[j])-Rbty(t[j]))^2)/L

                       }

errmse100 = klasma2*sum(XX)
errmse100

#-------------------------------------------
                  #Jeffreys#
#-------------------------------------------

Rb = 0

RJEF100  = function(x) {

for(i in 1:L) {

Rb = ((-T100+2*(1-exp(x)))/(-T100+(1-exp(x))))^n

              }
return(Rb)			   
 			     }


#-----------------
    #IMSE#
#-----------------

klasma1 = 1/L

klasma2 = 1/length(t)

XX = rep(NA,5)

for(j in 1:length(t)) {

XX[j] = sum((RJEF100(t[j])-Rbty(t[j]))^2)/L

                       }

errjef100 = klasma2*sum(XX)
errjef100

#-------------------------------------------
                  #Gamma#
#-------------------------------------------

-------------------------------------------
                 #gama=0.8#
#-------------------------------------------
gama=0.8
#------------------
#delta=0.5#
#------------------

delta=0.5

Rb = 0

RG1D1100  = function(x) {

for(i in 1:L) {

Rb = ((delta-T100+2*(1-exp(x)))/(delta-T100+(1-exp(x))))^(n+gama)

              }
return(Rb)			   
 			     }


#-----------------
    #IMSE#
#-----------------

klasma1 = 1/L

klasma2 = 1/length(t)

XX = rep(NA,5)

for(j in 1:length(t)) {

XX[j] = sum((RG1D1100(t[j])-Rbty(t[j]))^2)/L

                       }

errg1d1100 = klasma2*sum(XX)
errg1d1100

#------------------
#delta=3.0#
#------------------

delta=3.0

Rb = 0

RG1D2100  = function(x) {

for(i in 1:L) {

Rb = ((delta-T100+2*(1-exp(x)))/(delta-T100+(1-exp(x))))^(n+gama)

              }
return(Rb)			   
 			     }


#-----------------
    #IMSE#
#-----------------

klasma1 = 1/L

klasma2 = 1/length(t)

XX = rep(NA,5)

for(j in 1:length(t)) {

XX[j] = sum((RG1D2100(t[j])-Rbty(t[j]))^2)/L

                       }

errg1d2100 = klasma2*sum(XX)
errg1d2100

#-------------------------------------------
                 #gama=3.0#
#--------------------------------------------
gama=3.0
#------------------
#delta=0.5#
#------------------

delta=0.5

Rb = 0

RG2D1100  = function(x) {

for(i in 1:L) {

Rb = ((delta-T100+2*(1-exp(x)))/(delta-T100+(1-exp(x))))^(n+gama)

              }
return(Rb)			   
 			     }


#-----------------
    #IMSE#
#-----------------

klasma1 = 1/L

klasma2 = 1/length(t)

XX = rep(NA,5)

for(j in 1:length(t)) {

XX[j] = sum((RG2D1100(t[j])-Rbty(t[j]))^2)/L

                       }

errg2d1100 = klasma2*sum(XX)
errg2d1100

#------------------
#delta=3.0#
#------------------

delta=3.0

Rb = 0

RG2D2100  = function(x) {

for(i in 1:L) {

Rb = ((delta-T100+2*(1-exp(x)))/(delta-T100+(1-exp(x))))^(n+gama)

              }
return(Rb)			   
 			     }


#-----------------
    #IMSE#
#-----------------

klasma1 = 1/L

klasma2 = 1/length(t)

XX = rep(NA,5)

for(j in 1:length(t)) {

XX[j] = sum((RG2D2100(t[j])-Rbty(t[j]))^2)/L

                       }

errg2d2100 = klasma2*sum(XX)
errg2d2100

tb100 = c(errmle100,errumv100,errmse100,errjef100,errg1d1100,errg1d2100,errg2d1100,errg2d2100)
tb100


#-----------------------------------------
               #TABLE3#
#-----------------------------------------

TABLE3 = rbind(tb15,tb50,tb100)
colnames(TABLE3) = c('MLE','UMVUE','MSE','JEFFREYS'
                   ,'GAMMA(gama=0.8,delta=0.5)','GAMMA(gama=0.8,delta=3.0)'
                   ,'GAMMA(gama=3.0,delta=0.5)','GAMMA(gama=3.0,delta=3.0)')
rownames(TABLE3) = c('n=15','n=50','n=100')
TABLE3

















