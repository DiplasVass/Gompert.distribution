.libPaths("/home/testuser/Rlib")
install.packages('VGAM')
library(VGAM)

## Introduction ##

## Gompertz distribution

## Probability Density Function of Gompertz Distribution
 
t = seq(0,10,0.05)

gompd = function(phi,c,x) {

g = phi*exp(c*x+(phi*(1-exp(c*t))))

return(g)
}

plot(t,gompd(0.002,1,t),main='Prabability Density Function Of Gompertz Distribution (phi=0.002)')
lines(t,gompd(0.002,1,t),col='blue')


## Cumulative Density Function of Gompertz Distribution

cumgomp = function(phi,x) {

g = 1 - exp(phi*(1-exp(t)))

return(g)
}

plot(t,cumgomp(0.002,t),main='Cumylative Function Of Gompertz Distribution (phi=0.002)')
lines(t,cumgomp(0.002,t),col='red')


## Reliability Density Function of Gompertz Distribution

rely = function(phi,x) {

g = 1 - cumgomp(phi,x) 

return(g)
				}

plot(t,rely(0.002,t),main='Reliability Function Of Gompertz Distribution (phi=0.002)')
lines(t,rely(0.002,t),col='red')

------------------------------------------------------------------------------------------
                               ##Statistical Estimation##                                        
------------------------------------------------------------------------------------------

##### Estimators #####

-------------------------------------------------------------------------------------------

## Non - Bayes Estimators of the Shape Parameter Maximum Likelihood Estimator (MLE) ##

# We foumd the parameter estimation of phi with the help of the MLE 

T = sum(1-exp(t))

phimle = - (length(t)/T)
phimle

# The MLE for reliability function for the estimated phi (phihat) 

rmle = exp(phimle*(1-exp(t)))
rmle

plot(t,rely(phimle,t),col="red")
lines(t,rmle,col="green")

-------------------------------------------------------------------------------------------

## Uniformly Minimum Variance Unbiased Estimator (UMVUE) ##

# We foumd the parameter estimation of phi with the help of the UMVUE

T = sum(1-exp(t))

phiumv = - ((length(t)-1)/T)
phiumv

# The UMVUE for reliability function for the estimated phi (phihat) 

rumv = exp(phiumv*(1-exp(t)))
rumv

plot(t,rely(phiumv,t),col="red")
lines(t,rumv,col="green")

-------------------------------------------------------------------------------------------

## Minimum Mean Squared Error Estimators Method (MinMSE) ##

# We foumd the parameter estimation of phi with the help of the UMVUE

T = sum(1-exp(t))

phimse = - ((length(t)-2)/T)
phimse

# The UMVUE for reliability function for the estimated phi (phihat) 

rmse = exp(phimse*(1-exp(t)))
rmse

plot(t,rely(phimse,t),col="red")
lines(t,rmse,col="green")

--------------------------------------------------------------------------------------------

--------------------------------------------------------------------------------------------
                                  ##Bayes Estimation## 
--------------------------------------------------------------------------------------------


GAM = function(gama,delta) {

q = dgamma(t,alpha=gama,beta=1/delta) 
return(q)
                            }


### Bayes Estimation under (SISELF) with Gamma Prior ###

# We foumd the parameter estimation of phi

phibg = (n+gama-2)/(delta-T)
phibg

# The SISELF for reliability function for the estimated phi

A1  = delta - T + 2*(1-exp(t))
A2  = delta - T +   (1-exp(t))

rbg = (A1/A2)^(n+gama)


### Bayes Estimation under (SISELF) with Jefreys Prior ###

# We foumd the parameter estimation of phi

phibj = (n-2)/(-T)
phibj

# The SISELF for reliability function for the estimated phi

J1  =  - T + 2*(1-exp(t))
J2  =  - T +   (1-exp(t))

rbg = (J1/J2)^(n)

------------------------------------------------------------------------------------------------
                                    ##Simulation Study##
------------------------------------------------------------------------------------------------

----------------------------------------------------------------------------------------------
                                       # phi = 0.5 #
----------------------------------------------------------------------------------------------

L =5000 # We want the process to be repeated 5000 times 

########################### n = 15 ################################

n=15

## Construct a sample of Gompertz(0.5,1) with n = 15 observations

dat15 = rgompertz(n=15, scale = 1, shape=0.5)

## Monte Carlo Simulations (L = 5000) for n=15

Dat15 = replicate(L,rgompertz(15, scale = 1, shape=0.5))

-------------------------------------------------------------------------------------------
                        ##Statistical Estimation for phi##                                        
-------------------------------------------------------------------------------------------

##### Estimators #####

-------------------------------------------------------------------------------------------

# We foumd the T estimator which help us for the estimation of phi

T15 = 0

for(i in 1:L) {  # Find the 5000 values of T estimator  

T15[i] = sum(1-exp(Dat15[,i]))

		   }

T15

hist(T15)

---------------------------------------------------------------------------------------------

## Non - Bayes Estimators of the Shape Parameter Maximum Likelihood Estimator (MLE) ##

phimle15 = 0

for(i in 1:L) { # Find the 5000 values of phimle estimator

phimle15[i] = - (n/T15[i])

              }

phimle15

x15 = hist(phimle15,nclass=25)
plot(x15$breaks,c(x15$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for MLE method with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

-------------------------------------------------------------------------------------------

## Uniformly Minimum Variance Unbiased Estimator (UMVUE) ##

# We foumd the parameter estimation of phi with the help of the UMVUE

phiumv15 = 0

for(i in 1:L) { # Find the 5000 values of phiumv estimator

phiumv15[i] = - ((n-1)/T15[i])

		   }
phiumv15

y15 = hist(phiumv15,nclass=25)
plot(y15$breaks,c(y15$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for UMVUE method with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

-------------------------------------------------------------------------------------------

## Minimum Mean Squared Error Estimators Method (MinMSE) ##

# We foumd the parameter estimation of phi with the help of the MSE

phimse15 = 0 

for(i in 1:L) { # Find the 5000 values of phimse estimator

phimse15[i] = - ((n-2)/T15[i]) 

		  }
phimse15

z15 = hist(phimse15,nclass=25)
plot(z15$breaks,c(z15$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for MSE method with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')



## empirical values versus real value

r15 = cbind(mean(phimle15),mean(phiumv15),mean(phimse15),0.5)
r15

########################### n = 50 ################################

n=50

## Construct a sample of Gompertz(0.5,1) with n = 50 observations

dat50 = rgompertz(n=50, scale = 1, shape=0.5)

## Monte Carlo Simulations (L = 5000) for n=50

Dat50 = replicate(L,rgompertz(50, scale = 1, shape=0.5))

------------------------------------------------------------------------------------------
                        ##Statistical Estimation for phi##                                        
------------------------------------------------------------------------------------------

##### Estimators #####

-------------------------------------------------------------------------------------------

# We foumd the T estimator which help us for the estimation of phi

T50 = 0

for(i in 1:L) {  # Find the 5000 values of T estimator  

T50[i] = sum(1-exp(Dat50[,i]))

		   }

T50

hist(T50)

---------------------------------------------------------------------------------------------

## Non - Bayes Estimators of the Shape Parameter Maximum Likelihood Estimator (MLE) ##

phimle50 = 0

for(i in 1:L) { # Find the 5000 values of phimle estimator

phimle50[i] = - (n/T50[i])

              }

phimle50

x50 = hist(phimle50,nclass=25)
plot(x50$breaks,c(x50$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for MLE method with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

-------------------------------------------------------------------------------------------

## Uniformly Minimum Variance Unbiased Estimator (UMVUE) ##

# We foumd the parameter estimation of phi with the help of the UMVUE

phiumv50 = 0

for(i in 1:L) { # Find the 5000 values of phiumv estimator

phiumv50[i] = - ((n-1)/T50[i])

		   }
phiumv50

y50 = hist(phiumv50,nclass=25)
plot(y50$breaks,c(y50$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for UMVUE method with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

-------------------------------------------------------------------------------------------

## Minimum Mean Squared Error Estimators Method (MinMSE) ##

# We foumd the parameter estimation of phi with the help of the MSE

phimse50 = 0 

for(i in 1:L) { # Find the 5000 values of phimse estimator

phimse50[i] = - ((n-2)/T50[i]) 

		  }
phimse50

z50 = hist(phimse50,nclass=25)
plot(z50$breaks,c(z50$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for MSE method with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')



## empirical values versus real value

r50 = cbind(mean(phimle50),mean(phiumv50),mean(phimse50),0.5)
r50

########################### n = 100 ################################

n=100

## Construct a sample of Gompertz(0.5,1) with n = 50 observations

dat100 = rgompertz(n=100, scale = 1, shape=0.5)

## Monte Carlo Simulations (L = 5000) for n=50

Dat100 = replicate(L,rgompertz(100, scale = 1, shape=0.5))

------------------------------------------------------------------------------------------
                        ##Statistical Estimation for phi##                                        
------------------------------------------------------------------------------------------

##### Estimators #####

-------------------------------------------------------------------------------------------

# We foumd the T estimator which help us for the estimation of phi

T100 = 0

for(i in 1:L) {  # Find the 5000 values of T estimator  

T100[i] = sum(1-exp(Dat100[,i]))

		   }

T100

hist(T100)

---------------------------------------------------------------------------------------------

## Non - Bayes Estimators of the Shape Parameter Maximum Likelihood Estimator (MLE) ##

phimle100 = 0

for(i in 1:L) { # Find the 5000 values of phimle estimator

phimle100[i] = - (n/T100[i])

              }

phimle100

x100 = hist(phimle100,nclass=25)
plot(x100$breaks,c(x100$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for MLE method with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

-------------------------------------------------------------------------------------------

## Uniformly Minimum Variance Unbiased Estimator (UMVUE) ##

# We foumd the parameter estimation of phi with the help of the UMVUE

phiumv100 = 0

for(i in 1:L) { # Find the 5000 values of phiumv estimator

phiumv100[i] = - ((n-1)/T100[i])

		   }
phiumv100

y100 = hist(phiumv100,nclass=25)
plot(y100$breaks,c(y100$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for UMVUE method with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

-------------------------------------------------------------------------------------------

## Minimum Mean Squared Error Estimators Method (MinMSE) ##

# We foumd the parameter estimation of phi with the help of the MSE

phimse100 = 0 

for(i in 1:L) { # Find the 5000 values of phimse estimator

phimse100[i] = - ((n-2)/T100[i]) 

		  }
phimse100

z100 = hist(phimse100,nclass=25)
plot(z100$breaks,c(z100$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for MSE method with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

## empirical values versus real value

r100 = cbind(mean(phimle100),mean(phiumv100),mean(phimse100),0.5)
r100

## MATRIX 1 ##

matrix1 = rbind(r15,r50,r100)
rownames(matrix1) = c('n=15','n=50','n=100')
colnames(matrix1) = c('MLE','UMVUE','MSE','REAL')
matrix1

-------------------------------------
        #MSE OF ESTIMATORS#
-------------------------------------

-------------------------------------
            #METHODS#
-------------------------------------

-------------------------------------

############ n = 15 #################

-------------------------------------

-------------------------------------
              #MLE#
-------------------------------------

a=0
for(i in 1:L) {
a[i] = ((phimle15[i]-0.5)^2)
emle15 = sum(a)/L
		   }   

emle15 ## error of mle


-------------------------------------
              #UMVUE#
-------------------------------------

a=0
for(i in 1:L) {
a[i] = ((phiumv15[i]-0.5)^2)
eumv15 = sum(a)/L
		   }   

eumv15 ## error of umvue

-------------------------------------
              #MSE#
-------------------------------------

a=0
for(i in 1:L) {
a[i] = ((phimse15[i]-0.5)^2)
emse15 = sum(a)/L
		   }   

emse15 ## error of mse

-------------------------------------

############ n = 50 #################

-------------------------------------

-------------------------------------
              #MLE#
-------------------------------------

a=0
for(i in 1:L) {
a[i] = ((phimle50[i]-0.5)^2)
emle50 = sum(a)/L
		   }   

emle50 ## error of mle


-------------------------------------
              #UMVUE#
-------------------------------------

a=0
for(i in 1:L) {
a[i] = ((phiumv50[i]-0.5)^2)
eumv50 = sum(a)/L
		   }   

eumv50 ## error of umvue

-------------------------------------
              #MSE#
-------------------------------------

a=0
for(i in 1:L) {
a[i] = ((phimse50[i]-0.5)^2)
emse50 = sum(a)/L
		   }   

emse50 ## error of mse

-------------------------------------

############ n = 100 ################

-------------------------------------

-------------------------------------
              #MLE#
-------------------------------------

a=0
for(i in 1:L) {
a[i] = ((phimle100[i]-0.5)^2)
emle100 = sum(a)/L
		   }   

emle100 ## error of mle


-------------------------------------
              #UMVUE#
-------------------------------------

a=0
for(i in 1:L) {
a[i] = ((phiumv100[i]-0.5)^2)
eumv100 = sum(a)/L
		   }   

eumv100 ## error of umvue

-------------------------------------
              #MSE#
-------------------------------------

a=0
for(i in 1:L) {
a[i] = ((phimse100[i]-0.5)^2)
emse100 = sum(a)/L
		   }   

emse100 ## error of mse

## MATRIX 2 ##

er15 = cbind(emle15,eumv15,emse15)
er50 = cbind(emle50,eumv50,emse50)
er100 = cbind(emle100,eumv100,emse100)

matrix2 = rbind(er15,er50,er100)
rownames(matrix2) = c('n=15','n=50','n=100')
colnames(matrix2) = c('MLE','UMVUE','MSE')
matrix2


-------------------------------------------------------------------------------------------
                        ##Bayesian Estimation for phi##                                        
-------------------------------------------------------------------------------------------

##### Estimators #####

-------------------------------------------------------------------------------------------

---------------------------------------------
                  #gama = 0.8#
--------------------------------------------- 

gama=0.8

----------------------------------------------
                  #delta=0.5#
----------------------------------------------

delta=0.5

----------------------------------------------

GAM = function(gama,delta) {

q = dgamma(t,alpha=gama,beta=1/delta) 
return(q)
                            }

----------------------------------------------
                    #n = 15#
----------------------------------------------

n=15

### Bayes Estimation under (SISELF) with Gamma Prior ###

# We foumd the parameter estimation of phi

phibg15 = 0

for(i in 1:L) {

phibg15[i] = (n+gama-2)/(delta-T15[i])

		   }
phibg15

w15 = hist(phibg15,nclass=25)
plot(w15$breaks,c(w15$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for Gamma prior with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

### Bayes Estimation under (SISELF) with Jefreys Prior ###
P15 = - T15
# We foumd the parameter estimation of phi

phibj15 = 0

for(i in 1:L) {

phibj15[i] = (n-2)/(P15[i])

               }

phibj15

h15 = hist(phibj15,nclass=25)
plot(h15$breaks,c(h15$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for Jeffreys prior with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

b15 = cbind(mean(phibg15),mean(phibj15),0.5)
b15

----------------------------------------------
                    #n = 50#
----------------------------------------------

n=50

### Bayes Estimation under (SISELF) with Gamma Prior ###

# We foumd the parameter estimation of phi

phibg50 = 0

for(i in 1:L) {

phibg50[i] = (n+gama-2)/(delta-T50[i])

		   }
phibg50

w50 = hist(phibg50,nclass=25)
plot(w50$breaks,c(w50$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for Gamma prior with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

### Bayes Estimation under (SISELF) with Jefreys Prior ###
P50 = - T50
# We foumd the parameter estimation of phi

phibj50 = 0

for(i in 1:L) {

phibj50[i] = (n-2)/(P50[i])

               }

phibj50

h50 = hist(phibj50,nclass=25)
plot(h50$breaks,c(h50$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for Jeffreys prior with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

b50 = cbind(mean(phibg50),mean(phibj50),0.5)
b50

----------------------------------------------
                    #n = 100#
----------------------------------------------

n=100

### Bayes Estimation under (SISELF) with Gamma Prior ###

# We foumd the parameter estimation of phi

phibg100 = 0

for(i in 1:L) {

phibg100[i] = (n+gama-2)/(delta-T100[i])

		   }
phibg100

w100 = hist(phibg100,nclass=25)
plot(w100$breaks,c(w100$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for Gamma prior with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

### Bayes Estimation under (SISELF) with Jefreys Prior ###
P100 = - T100
# We foumd the parameter estimation of phi

phibj100 = 0

for(i in 1:L) {

phibj100[i] = (n-2)/(P100[i])

               }

phibj100

h100 = hist(phibj100,nclass=25)
plot(h100$breaks,c(h100$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for Jeffreys prior with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

b100 = cbind(mean(phibg100),mean(phibj100),0.5)
b100

## MATRIX 3i ##

matrix3i = rbind(b15,b50,b100)
rownames(matrix3i) = c('n=15','n=50','n=100')
colnames(matrix3i) = c('Gamma','Jeffreys','REAL')
matrix3i

-------------------------------------
        #MSE OF ESTIMATORS#
-------------------------------------

-------------------------------------
            #METHODS#
-------------------------------------

-------------------------------------

############ n = 15 #################

-------------------------------------

-------------------------------------
           #Gamma Prior#
-------------------------------------

a=0
for(i in 1:L) {
a[i] = ((phibg15[i]-0.5)^2)
ebg15 = sum(a)/L
		   }   

ebg15 ## error of gamma prior

-------------------------------------

############ n = 50 #################

-------------------------------------

-------------------------------------
           #Gamma Prior#
-------------------------------------

a=0
for(i in 1:L) {
a[i] = ((phibg50[i]-0.5)^2)
ebg50 = sum(a)/L
		   }   

ebg50 ## error of gamma prior

-------------------------------------

############ n = 100 ################

-------------------------------------

-------------------------------------
           #Gamma Prior#
-------------------------------------

a=0
for(i in 1:L) {
a[i] = ((phibg100[i]-0.5)^2)
ebg100 = sum(a)/L
		   }   

ebg100 ## error of gamma prior

## MATRIX 3ii ##

matrix3ii = rbind(ebg15,ebg50,ebg100)
rownames(matrix3ii) = c('n=15','n=50','n=100')
colnames(matrix3ii) = c('Gamma Prior (delta = 0.5)')
matrix3ii

----------------------------------------------
                  #delta=3.0#
----------------------------------------------

delta=3

----------------------------------------------

GAM = function(gama,delta) {

q = dgamma(t,alpha=gama,beta=1/delta) 
return(q)
                            }

----------------------------------------------
                    #n = 15#
----------------------------------------------

n=15

### Bayes Estimation under (SISELF) with Gamma Prior ###

# We foumd the parameter estimation of phi

phibg15 = 0

for(i in 1:L) {

phibg15[i] = (n+gama-2)/(delta-T15[i])

		   }
phibg15

w15 = hist(phibg15,nclass=25)
plot(w15$breaks,c(w15$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for Gamma prior with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

### Bayes Estimation under (SISELF) with Jefreys Prior ###
P15 = - T15
# We foumd the parameter estimation of phi

phibj15 = 0

for(i in 1:L) {

phibj15[i] = (n-2)/(P15[i])

               }

phibj15

h15 = hist(phibj15,nclass=25)
plot(h15$breaks,c(h15$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for Jeffreys prior with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

b15 = cbind(mean(phibg15),mean(phibj15),0.5)
b15

----------------------------------------------
                    #n = 50#
----------------------------------------------

n=50

### Bayes Estimation under (SISELF) with Gamma Prior ###

# We foumd the parameter estimation of phi

phibg50 = 0

for(i in 1:L) {

phibg50[i] = (n+gama-2)/(delta-T50[i])

		   }
phibg50

w50 = hist(phibg50,nclass=25)
plot(w50$breaks,c(w50$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for Gamma prior with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

### Bayes Estimation under (SISELF) with Jefreys Prior ###
P50 = - T50
# We foumd the parameter estimation of phi

phibj50 = 0

for(i in 1:L) {

phibj50[i] = (n-2)/(P50[i])

               }

phibj50

h50 = hist(phibj50,nclass=25)
plot(h50$breaks,c(h50$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for Jeffreys prior with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

b50 = cbind(mean(phibg50),mean(phibj50),0.5)
b50

----------------------------------------------
                    #n = 100#
----------------------------------------------

n=100

### Bayes Estimation under (SISELF) with Gamma Prior ###

# We foumd the parameter estimation of phi

phibg100 = 0

for(i in 1:L) {

phibg100[i] = (n+gama-2)/(delta-T100[i])

		   }
phibg100

w100 = hist(phibg100,nclass=25)
plot(w100$breaks,c(w100$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for Gamma prior with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

### Bayes Estimation under (SISELF) with Jefreys Prior ###
P100 = - T100
# We foumd the parameter estimation of phi

phibj100 = 0

for(i in 1:L) {

phibj100[i] = (n-2)/(P100[i])

               }

phibj100

h100 = hist(phibj100,nclass=25)
plot(h100$breaks,c(h100$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for Jeffreys prior with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

b100 = cbind(mean(phibg100),mean(phibj100),0.5)
b100

## MATRIX 3iii ##

matrix3iii = rbind(b15,b50,b100)
rownames(matrix3iii) = c('n=15','n=50','n=100')
colnames(matrix3iii) = c('Gamma','Jeffreys','REAL')
matrix3iii

-------------------------------------
        #MSE OF ESTIMATORS#
-------------------------------------

-------------------------------------
            #METHODS#
-------------------------------------

-------------------------------------

############ n = 15 #################

-------------------------------------

-------------------------------------
           #Gamma Prior#
-------------------------------------

a=0
for(i in 1:L) {
a[i] = ((phibg15[i]-0.5)^2)
ebg15 = sum(a)/L
		   }   

ebg15 ## error of gamma prior

-------------------------------------

############ n = 50 #################

-------------------------------------

-------------------------------------
           #Gamma Prior#
-------------------------------------

a=0
for(i in 1:L) {
a[i] = ((phibg50[i]-0.5)^2)
ebg50 = sum(a)/L
		   }   

ebg50 ## error of gamma prior

-------------------------------------

############ n = 100 ################

-------------------------------------

-------------------------------------
           #Gamma Prior#
-------------------------------------

a=0
for(i in 1:L) {
a[i] = ((phibg100[i]-0.5)^2)
ebg100 = sum(a)/L
		   }   

ebg100 ## error of gamma prior

## MATRIX 3iv ##

matrix3iv = rbind(ebg15,ebg50,ebg100)
rownames(matrix3iv) = c('n=15','n=50','n=100')
colnames(matrix3iv) = c('Gamma Prior (delta = 3.0)')
matrix3iv

---------------------------------------------
                  #gama = 3.0#
--------------------------------------------- 

gama=0.8

----------------------------------------------
                  #delta=0.5#
----------------------------------------------

delta=0.5

----------------------------------------------

GAM = function(gama,delta) {

q = dgamma(t,alpha=gama,beta=1/delta) 
return(q)
                            }

----------------------------------------------
                    #n = 15#
----------------------------------------------

n=15

### Bayes Estimation under (SISELF) with Gamma Prior ###

# We foumd the parameter estimation of phi

phibg15 = 0

for(i in 1:L) {

phibg15[i] = (n+gama-2)/(delta-T15[i])

		   }
phibg15

w15 = hist(phibg15,nclass=25)
plot(w15$breaks,c(w15$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for Gamma prior with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

### Bayes Estimation under (SISELF) with Jefreys Prior ###
P15 = - T15
# We foumd the parameter estimation of phi

phibj15 = 0

for(i in 1:L) {

phibj15[i] = (n-2)/(P15[i])

               }

phibj15

h15 = hist(phibj15,nclass=25)
plot(h15$breaks,c(h15$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for Jeffreys prior with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

b15 = cbind(mean(phibg15),mean(phibj15),0.5)
b15

----------------------------------------------
                    #n = 50#
----------------------------------------------

n=50

### Bayes Estimation under (SISELF) with Gamma Prior ###

# We foumd the parameter estimation of phi

phibg50 = 0

for(i in 1:L) {

phibg50[i] = (n+gama-2)/(delta-T50[i])

		   }
phibg50

w50 = hist(phibg50,nclass=25)
plot(w50$breaks,c(w50$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for Gamma prior with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

### Bayes Estimation under (SISELF) with Jefreys Prior ###
P50 = - T50
# We foumd the parameter estimation of phi

phibj50 = 0

for(i in 1:L) {

phibj50[i] = (n-2)/(P50[i])

               }

phibj50

h50 = hist(phibj50,nclass=25)
plot(h50$breaks,c(h50$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for Jeffreys prior with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

b50 = cbind(mean(phibg50),mean(phibj50),0.5)
b50

----------------------------------------------
                    #n = 100#
----------------------------------------------

n=100

### Bayes Estimation under (SISELF) with Gamma Prior ###

# We foumd the parameter estimation of phi

phibg100 = 0

for(i in 1:L) {

phibg100[i] = (n+gama-2)/(delta-T100[i])

		   }
phibg100

w100 = hist(phibg100,nclass=25)
plot(w100$breaks,c(w100$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for Gamma prior with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

### Bayes Estimation under (SISELF) with Jefreys Prior ###
P100 = - T100
# We foumd the parameter estimation of phi

phibj100 = 0

for(i in 1:L) {

phibj100[i] = (n-2)/(P100[i])

               }

phibj100

h100 = hist(phibj100,nclass=25)
plot(h100$breaks,c(h100$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for Jeffreys prior with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

b100 = cbind(mean(phibg100),mean(phibj100),0.5)
b100

## MATRIX 4i ##

matrix4i = rbind(b15,b50,b100)
rownames(matrix4i) = c('n=15','n=50','n=100')
colnames(matrix4i) = c('Gamma','Jeffreys','REAL')
matrix4i

-------------------------------------
        #MSE OF ESTIMATORS#
-------------------------------------

-------------------------------------
            #METHODS#
-------------------------------------

-------------------------------------

############ n = 15 #################

-------------------------------------

-------------------------------------
           #Gamma Prior#
-------------------------------------

a=0
for(i in 1:L) {
a[i] = ((phibg15[i]-0.5)^2)
ebg15 = sum(a)/L
		   }   

ebg15 ## error of gamma prior

-------------------------------------

############ n = 50 #################

-------------------------------------

-------------------------------------
           #Gamma Prior#
-------------------------------------

a=0
for(i in 1:L) {
a[i] = ((phibg50[i]-0.5)^2)
ebg50 = sum(a)/L
		   }   

ebg50 ## error of gamma prior

-------------------------------------

############ n = 100 ################

-------------------------------------

-------------------------------------
           #Gamma Prior#
-------------------------------------

a=0
for(i in 1:L) {
a[i] = ((phibg100[i]-0.5)^2)
ebg100 = sum(a)/L
		   }   

ebg100 ## error of gamma prior

## MATRIX 4ii ##

matrix4ii = rbind(ebg15,ebg50,ebg100)
rownames(matrix4ii) = c('n=15','n=50','n=100')
colnames(matrix4ii) = c('Gamma Prior (delta = 0.5)')
matrix4ii

----------------------------------------------
                  #delta=3.0#
----------------------------------------------

delta=3

----------------------------------------------

GAM = function(gama,delta) {

q = dgamma(t,alpha=gama,beta=1/delta) 
return(q)
                            }

----------------------------------------------
                    #n = 15#
----------------------------------------------

n=15

### Bayes Estimation under (SISELF) with Gamma Prior ###

# We foumd the parameter estimation of phi

phibg15 = 0

for(i in 1:L) {

phibg15[i] = (n+gama-2)/(delta-T15[i])

		   }
phibg15

w15 = hist(phibg15,nclass=25)
plot(w15$breaks,c(w15$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for Gamma prior with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

### Bayes Estimation under (SISELF) with Jefreys Prior ###
P15 = - T15
# We foumd the parameter estimation of phi

phibj15 = 0

for(i in 1:L) {

phibj15[i] = (n-2)/(P15[i])

               }

phibj15

h15 = hist(phibj15,nclass=25)
plot(h15$breaks,c(h15$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for Jeffreys prior with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

b15 = cbind(mean(phibg15),mean(phibj15),0.5)
b15

----------------------------------------------
                    #n = 50#
----------------------------------------------

n=50

### Bayes Estimation under (SISELF) with Gamma Prior ###

# We foumd the parameter estimation of phi

phibg50 = 0

for(i in 1:L) {

phibg50[i] = (n+gama-2)/(delta-T50[i])

		   }
phibg50

w50 = hist(phibg50,nclass=25)
plot(w50$breaks,c(w50$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for Gamma prior with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

### Bayes Estimation under (SISELF) with Jefreys Prior ###
P50 = - T50
# We foumd the parameter estimation of phi

phibj50 = 0

for(i in 1:L) {

phibj50[i] = (n-2)/(P50[i])

               }

phibj50

h50 = hist(phibj50,nclass=25)
plot(h50$breaks,c(h50$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for Jeffreys prior with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

b50 = cbind(mean(phibg50),mean(phibj50),0.5)
b50

----------------------------------------------
                    #n = 100#
----------------------------------------------

n=100

### Bayes Estimation under (SISELF) with Gamma Prior ###

# We foumd the parameter estimation of phi

phibg100 = 0

for(i in 1:L) {

phibg100[i] = (n+gama-2)/(delta-T100[i])

		   }
phibg100

w100 = hist(phibg100,nclass=25)
plot(w100$breaks,c(w100$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for Gamma prior with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

### Bayes Estimation under (SISELF) with Jefreys Prior ###
P100 = - T100
# We foumd the parameter estimation of phi

phibj100 = 0

for(i in 1:L) {

phibj100[i] = (n-2)/(P100[i])

               }

phibj100

h100 = hist(phibj100,nclass=25)
plot(h100$breaks,c(h100$counts,0),type="s",xlab="Data",ylab="Frequency",
     main='Histogram of phi values for Jeffreys prior with 
         the real phi parameter as blue line')
abline(v=0.5,lty=2,col='blue')

b100 = cbind(mean(phibg100),mean(phibj100),0.5)
b100

## MATRIX 4iii ##

matrix4iii = rbind(b15,b50,b100)
rownames(matrix4iii) = c('n=15','n=50','n=100')
colnames(matrix4iii) = c('Gamma','Jeffreys','REAL')
matrix4iii

-------------------------------------
        #MSE OF ESTIMATORS#
-------------------------------------

-------------------------------------
            #METHODS#
-------------------------------------

-------------------------------------

############ n = 15 #################

-------------------------------------

-------------------------------------
           #Gamma Prior#
-------------------------------------

a=0
for(i in 1:L) {
a[i] = ((phibg15[i]-0.5)^2)
ebg15 = sum(a)/L
		   }   

ebg15 ## error of gamma prior

-------------------------------------

############ n = 50 #################

-------------------------------------

-------------------------------------
           #Gamma Prior#
-------------------------------------

a=0
for(i in 1:L) {
a[i] = ((phibg50[i]-0.5)^2)
ebg50 = sum(a)/L
		   }   

ebg50 ## error of gamma prior

-------------------------------------

############ n = 100 ################

-------------------------------------

-------------------------------------
           #Gamma Prior#
-------------------------------------

a=0
for(i in 1:L) {
a[i] = ((phibg100[i]-0.5)^2)
ebg100 = sum(a)/L
		   }   

ebg100 ## error of gamma prior

## MATRIX 4iv ##

matrix4iv = rbind(ebg15,ebg50,ebg100)
rownames(matrix4iv) = c('n=15','n=50','n=100')
colnames(matrix4iv) = c('Gamma Prior (delta = 3.0)')
matrix4iv

############ n = 15 #################

-------------------------------------

-------------------------------------
           #Jeffreys Prior#
-------------------------------------

a=0
for(i in 1:L) {
a[i] = ((phibj15[i]-0.5)^2)
ebj15 = sum(a)/L
		   }   

ebj15 ## error of jeffreys prior

-------------------------------------

############ n = 50 #################

-------------------------------------

-------------------------------------
           #Jeffreys Prior#
-------------------------------------

a=0
for(i in 1:L) {
a[i] = ((phibj50[i]-0.5)^2)
ebj50 = sum(a)/L
		   }   

ebj50 ## error of jeffreys prior

-------------------------------------

############ n = 100 ################

-------------------------------------

-------------------------------------
           #Jeffreys Prior#
-------------------------------------

a=0
for(i in 1:L) {
a[i] = ((phibj100[i]-0.5)^2)
ebj100 = sum(a)/L
		   }   

ebj100 ## error of jeffreys prior

## MATRIX 5 ##

matrix5 = rbind(ebj15,ebj50,ebj100)
rownames(matrix5) = c('n=15','n=50','n=100')
colnames(matrix5) = c('Jeffreys Prior')
matrix5
matrix3i[,2]
#TABLE 1#

#co1.i : (i=1,2,3) values of phi with mle,umvue,mse (expected and mse)

co1.1 = rbind(matrix1[1,1],matrix2[1,1]
           ,matrix1[2,1],matrix2[2,1]
           ,matrix1[3,1],matrix2[3,1])

co1.2 = rbind(matrix1[1,2],matrix2[1,2]
           ,matrix1[2,2],matrix2[2,2]
           ,matrix1[3,2],matrix2[3,2])

co1.3 = rbind(matrix1[1,3],matrix2[1,3]
           ,matrix1[2,3],matrix2[2,3]
           ,matrix1[3,3],matrix2[3,3])

#co2 : values of phi with Jeffreys prior
co2 = rbind(matrix3iii[1,2],matrix5[1]
           ,matrix3iii[2,2],matrix5[2]
           ,matrix3iii[3,2],matrix5[3]) 

#co3 : values of phi with Gamma prior (gama=0.8,delta=0.5)
co3 = rbind(matrix3i[1,1],matrix3ii[1,1]
           ,matrix3i[2,1],matrix3ii[2,1]
           ,matrix3i[3,1],matrix3ii[3,1])

#co4 : values of phi with Gamma prior (gama=0.8,delta=3.0)
co4 = rbind(matrix3iii[1,1],matrix3iv[1,1]
           ,matrix3iii[2,1],matrix3iv[2,1]
           ,matrix3iii[3,1],matrix3iv[3,1])

#co5 : values of phi with Gamma prior (gama=3.0,delta=0.5)
co5 = rbind(matrix4i[1,1],matrix4ii[1,1]
           ,matrix4i[2,1],matrix4ii[2,1]
           ,matrix4i[3,1],matrix4ii[3,1])

#co6 : values of phi with Gamma prior (gama=3.0,delta=3.0)
co6 = rbind(matrix4iii[1,1],matrix4iv[1,1]
           ,matrix4iii[2,1],matrix4iv[2,1]
           ,matrix4iii[3,1],matrix4iv[3,1])

cret = c('EXP','MSE','EXP','MSE','EXP','MSE')



TABLE1 = matrix(c(co1.1,co1.2,co1.3,co2,co3,co4,co5,co5,co6),6)
rownames(TABLE1) = c('n=15  EXP','n=15  MSE','n=50  EXP','n=50  MSE','n=100 EXP','n=100 MSE') 
colnames(TABLE1) = c('n   Criteria','MLE','UMVUE','MSE','Jeffreys prior'
                     ,'Gamma prior (gama=0.8,delta=0.5)','Gamma prior (gama=0.8,delta=3.0)'
                     ,'Gamma prior (gama=3.0,delta=0.5)','Gamma prior (gama=3.0,delta=3.0)')

TABLE1



























 
