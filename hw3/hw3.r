'''
In this homework, you are supposed to learn how to simulate and estimate AR/MA models
in R, and how to use ACF and PACF to identify models
'''

# 3.1
'''
Consider an AR(1) process f(t) = ø*f(t-1) + z(t), z~N(0,4)   |ø| < 1
Simulate with ø = 0.7, t=1..n, n=10000=10^4

Note: ACF of AR(1): y = ø^h, h >= 1
partial autocorrelation function alpha(h) = 1 if h=1 else 0
'''

#A: Plot the sample autocorrelation function ρ_hat sample partial autocorrelation function α_hat of the simulated process.
#How are these consistent with an AR(1) process? Can you see the relation between ø and ρ_hat(1)


set.seed(1234)
y <- arima.sim(
	model = list(ar = 0.7),
	n = 10000,
	innov = rnorm(10000, sd = 2)
)

par(mfrow=c(1,2)) #create a window with 1 row and 2 columns
AutoCorrelation = acf(x, plot = FALSE)
plot(AutoCorrelation, main = "Auto Corr. Func. AR(1)")
PartialAutoCorr = pacf(x, plot=FALSE)
plot(PartialAutoCorr, main="Partial Auto Corr. Func.")

#The ACF "tails off" to || x_n - x_m || < Epsilon after t>=10
#The PACF "cuts off" at t=2

'''
f\Mdl|      AR(P)     | MA(q)           | ARMA(p,q)
----------------------------------------------
ACF  |Tails off       | cuts off at t=q | Tails off
PACF |cuts off at t=p | Tails off       | Tails off
'''

alpha_1 =pacf(x, plot=FALSE)$acf[1]

#[1] 0.7011247
# alpha(1) approx ø

#B: Make a plot comparing the sample autocorrelation with the theoretical autocorrelation. Do the same for the partial autocorrelation

phi = .7
acf_theo = phi^(0:40) # Numerical Theoretical ACF with up to 40 lags - length 41 y(h) = ø^h , h >= 1
pacf_theo = c(phi, rep(0,39)) # Theoretical PACF with up to 40 lags - length 40: p(h) = ø^h , if h=1 else p(h)=0
# Plotting ACF and PACF
par(mfrow=c(1,2))

acf_empirical = acf(x, plot=FALSE)$acf
pacf_empirical = pacf(x, plot=FALSE)$acf

# length acf(x, plot=FALSE)$acf = 41
plot(0:40, acf_empirical, type = "h", ylab = "ACF", xlab = "h")
lines(0:40, acf_theo , col=2, type = "l", lty=2)
legend("topright", col = 1:2, legend = c(expression(hat(rho)(h)),expression(rho(h))), lty=1:2)

# length( pacf(x,plot=FALSE)$acf ) = 40
plot(1:40, pacf_empirical, type = "h", ylab = "PACF", xlab ="h")
lines(1:40, pacf_theo, col=2, type = "l",lty=2)
legend("topright", col = 1:2, legend = c(expression(hat(alpha)(h)),expression(alpha(h))), lty=1:2)

# length [acf = pacf +1]

# The theoretical and empirical seem to agree.

#C: Estimate an AR(1) and an AR(2) model for the simulated data using the arima function in R. 
#   Compare the two models. 
#   Is ø_hat_2 signficantly different from zero? 
#   Compare the quality of ø_1 for the two models in terms of bias and standard error.

arima(x, order =c(1,0,0), include.mean = FALSE) # AR(1) # order (p,d,q)(P,Q,D) - (non-seasonal) component

'''
Call:
arima(x = x, order = c(1, 0, 0), include.mean = FALSE)

Coefficients:
         ar1
      0.7012 = Ø
s.e.  0.0071 = standard error

sigma^2 estimated as 3.901:  log likelihood = -20996.39,  aic = 41996.78
'''

arima(x, order =c(2,0,0), include.mean = FALSE) # AR(2)

'''
Call:
arima(x = x, order = c(2, 0, 0), include.mean = FALSE)

Coefficients:
         ar1      ar2
      0.7023  -0.0015
s.e.  0.0100   0.0100

sigma^2 estimated as 3.901:  log likelihood = -20996.38,  aic = 41998.76
'''

# Q:Is ø_hat_2 signficantly different from zero? 
# A: ø_hat_2 = -0.0015 - close to zero

# Q: Compare the quality of ø_hat_1 for the two models in terms of bias and standard error.
# A: Coefficients for AR(1) -> ø_hat_1=0.712, and AR(2) -> ø_hat_1 =0.7023 There is a small difference
#    Standard error - AR(1) -> s.e    =0.0071, and AR(2) ->s.e.    =0.01   There is a small difference

#Note: ø_hat_2 =0.01 and s.e for it is 0.01 so ø_hat_2 is close to zero
#
# AR(1) =>  log likelihood = -20996.39
# AR(2) =>  log likelihood = -20996.38

# 3.2 
# C:  Simulate the  process as an AR(1) model 
#     with θ = 0.7 for t = 1, . . . , n with n = 10000 by implementing MA(1)
# Y_t = θ Z_t−1 + Z

y = arima.sim(
	model = list(ma = c(.7)),
	n = 10000,
	innov = rnorm(10000, sd = 2)
)

#D: Plot the sample ACF and PACF for the simulated data. How are these consistent with an MA(1) process
par(mfrow=c(1,2))
acf(y)
pacf(y)

#ACF cuts off after one lag 
#PACF tails off - like MA(p) model

# E: Fit an MA(a) and an AR(5) model to the simulated data. 
#    Compare the respective AR coefficients with the correspodning powers of θ in (2). 
#    How well does it fit?
#    Will an AR(10) improve the result?

theta = 0.7
ma1 = arima(y, order = c(0,0,1), include.mean=FALSE) #MA(1)
'''
Call:
arima(x = y, order = c(0, 0, 1), include.mean = FALSE)

Coefficients:
         ma1
      0.6987
s.e.  0.0072

sigma^2 estimated as 3.992:  log likelihood = -21111.36,  aic = 42226.73
'''

ar5 = arima(y, order = c(5,0,0), include.mean=FALSE)
'''
Call:
arima(x = y, order = c(5, 0, 0), include.mean = FALSE)

Coefficients:
         ar1      ar2     ar3      ar4     ar5
      0.6876  -0.4605  0.2956  -0.1814  0.0764
s.e.  0.0100   0.0120  0.0125   0.0120  0.0100

sigma^2 estimated as 4.028:  log likelihood = -21155.74,  aic = 42323.49
'''

# From (2) we should have:
'''
> -0.7^2
[1] -0.49
> 0.7^3
[1] 0.343
> -0.7^4
[1] -0.2401
> 0.7^5
[1] 0.16807
'''
# For j=1..4 the terms match +-0.05 but is more for j=5.

ar10 = arima(y, order = c(10,0,0), include.mean=FALSE)
rbind(ar10$coef,theta^(1:10)*(-1)^(2:11))

#The terms match very well

