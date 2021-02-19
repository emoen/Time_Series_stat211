# 2.2 a)
'''
Install the R-package astsa, load the package and the dataset varve. You may use
the following code:
	install.packages("astsa")
	library(astsa)
	data(varve)
Plot the glacial varve data.
The time series exhibits some nonstationarity that can be improved by transforming to 
logarithms and some additional nonstationarity that can be corrected by
differencing the logarithms.
'''

install.packages("astsa")
library(astsa)
data(varve)

plot(varve)

#b)
'''
Argue that the glacial varves series, say Xt, 
exhibits heteroscedasticity by computing the sample variance over the first half and the second half of the data. 
Argue that the transformation Yt = log Xt stabilizes the variance over the series. 
Plot the histograms of Xt and Yt to see whether the approximation to normality 
is improved by transforming the data.
'''

hist(varve, breaks = 10)
#Not so different
for (i in c(0:9)){
	sample = sample.int(n=length(varve), size=0.5*length(varve))
	print(var(varve[sample]))
	print(var(varve[-sample]))
	print("********")
}

#different variance
subset = seq(from = 1,to = length(varve)/2, by = 1)
var(varve[subset])
var(varve[-subset])

#log-transformation
y = log(varve)
par(mfrow=c(1,2))
hist(y, breaks = 20, prob=TRUE)
lines(dt <- seq(0,150,0.1), dnorm(dt, mean(y), sd(y)), col=2, lty=2)
hist(varve, breaks=10, prob=TRUE)
lines(dt <- seq(0,150,0.1), dnorm(dt, mean(varve),sd(varve)), col=2, lty=2)

'''
heteroscedasticity is the variability of a variable is unequal across the range of values of a second variable that predicts it
After log transform - the variable is symmetic about its mean - so more normal.
'''

#c) Plot the series Yt
plot(log(varve),lty=1)

#d) Examine the sample ACF of Yt and comment.
acf(log(varve), lag.max= 40))
acf(log(varve), lag.max= 150))
acf(varve, lag.max= 150))

#looking at plot in c) Y_t oscilates, acf(y) oscilates and acf is declining 
#So there is trend and seasonality

#e)
'''
Compute the difference Ut = Yt - Yt−1 , examine its time plot and sample ACF,
and argue that differencing the logged varve data produces a reasonably stationary
series. Can you think of a practical interpretation for Ut ?
'''

#U_t is the differenced (differenciated) time-series of Y_t - so without trend
U_t_x= diff(varve, lag = 1)
plot(U_t_x)
U_t= diff(y, lag = 1)
plot(U_t)

acf(U_t, lag.max=150)

#U_t is the definition of derivative when h=1. Since time discrete this is smallest h.

#f)
'''
Based on the sample ACF of the differenced transformed series computed in (d),
argue that the model in (1) below might be reasonable. 

Assume U_t = mu + Z_t + theta*Z_t-1

is stationary when the inputs Zt are assumed independent with mean 0 and variance sigma_z**2
Then, show ACVF: (1+theta**2) , h=0, theta*sima_z**2, h=1, 0 else
'''

#see paper solution - its a MA(1) processs

# g)
'''
Based on part (f), use ρbU (1) and the estimate of the variance of Ut
, γ_U (1) ), to
derive estimates of θ and σ
2
Z
. This is an application of the method of moments
from classical statistics, where estimators of the parameters are derived by equating
sample moments to theoretical moments. You can calculate the empirical autocorrelation and autocovariance functions 
using the acf function, by using the following syntax
	acf(u, type = "correlation")
	acf(u, type = "covariance")
'''

acf(U_t, type = "correlation")
acf(U_t, type = "covariance")

(a<-acf(U_t, plot=FALSE)$acf[2])
(b<-acf(U_t, type="covariance", plot=FALSE)$acf[1])
(theta <- (1 + c(-1,1) * sqrt(1 - 4 * a^2))/(2 * a))
(sigma2<- b / (1 + theta^2))
arima(U_t, order=c(0,0,1), include.mean=FALSE)

########
p <- as.vector(acf(Ut, type = "correlation", lag.max = 1, plot=FALSE)$acf[2]) #rho
g <- as.vector(acf(Ut, type = "covariance", lag.max = 0, plot=FALSE)$acf) #g=gamma

theta <- c((1+sqrt(1-4*p^2))/(2*p),(1-sqrt(1-4*p^2))/(2*p))
sigmasq <- g/(1+theta^2)

theta
sigmasq
#Need |theta| < 1 to be invertible, thus the solution is
#theta hat = -0.494 and sigmasq hat = 0.266



