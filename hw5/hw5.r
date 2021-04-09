# 5.1
# Let q = 5. Draw θj
# from the uniform distribution on [−2, 2] for j = 1, . . . , q. 
# Let σ^2_z= [1 + sum(θ^2_j)]^-1 and 
# {Zt} dist WN(0, σ2Z) and define the MA(q) process:
# a) Calculate and plot the autocovariance function defined by {(θj, j = 1, . . . , q), σ2Z} .

set.seed(1234)
q = 5
#draw U[-2,2] q times
theta = runif(n=q, min = -2, max = 2);theta
## [1]  -1.5451864  0.4891976  0.4370989  0.4935178  1.4436615
sigma = sqrt( 1/(1+sum(theta^2)) );sigma
## [1] 0.4033803


# Auto-covariance function for MA(q)
# %*% is matrix multiplication and sum
gammaX = function(h, theta, sigma) {
    h = abs(h)
    if ( h > length(theta)) return(0)	
	theta_k = c(1,theta)
    theta_kh = c(theta_k[(h+1):length(theta_k)], rep(0,h ))
    g = sigma^2 * sum(theta_k * theta_kh ) #The formula for ACVF of MA(q)
    return (g)
}


h = 0:10
'''
sapply(X, FUN)
Arguments:
-X: A vector or an object
-FUN: Function applied to each element of x
'''
gam = sapply(h, gammaX, theta=theta, sigma=sigma)
print(gam)
plot(h, gam, type="h", col=4)

#aacvf 
library("itsmr")
a = specify(ar=c(0),ma=theta, sigma2=sigma^2)
gamma.itsmr = aacvf(a,10)

plot(h, gamma.itsmr, type="h")
gamma.itsmr


# b) Make an R program for IA and find (θ∞, ν∞) = {(θ∞j, j = 1, . . . , q), ν∞}
N = 100
q = 5
IA = function(N, gamma.vec){
	#I) Initialisation page 30
	v = numeric(N) #initialize nu
	gamma = c(gamma.vec[1:(q+1)], rep(0,N)) #initialize
	Theta = matrix(0, ncol=N, nrow=N) #initialize
	v[1] = gamma[1] #gamma0
	Theta[1,1] <- gamma[2]/v[1] #theta1,1=gamma(1)/nu0
	v[2] = (1-Theta[1,1]^2)*v[1] #nu1 = (1-theta1,1^2)nu0
	
	#II) Start up part, for n=2,..,q
	for(n in 2:q) {
		Theta[n,n] = gamma[n+1]/v[1] #theta(n,n) = gamma(n)/nu0
		for (k in 1:(n-1)) {
			Theta[n,n-k] = (gamma[n+1-k]-sum(Theta[k, k:1]*Theta[n, n:(n+1-k)]*v[1:k]))/v[k+1] 
		}
		v[n+1] = gamma[1] - sum( Theta[n,n:1]^2*v[1:n] ) #nu(n)
	}	

	#III) Steady state part; for n>=q+1, formula on IA slide page 31
	for(n in (q+1):N) {
		Theta[n,q] = gamma[q+1]/v[n+1-q] #Theta(n,q)=gamma(q)/nu(n-q)
		for(k in 1:q-1){
			Theta[n,q-k] = (gamma[q+1-k] - sum(Theta[n-q+k,k:1]*Theta[n,q:(q-k+1)]*v[(n-q):(n-q+k-1)]))/v[n-q+k]
		}
		v[n+1] = gamma[1] - sum(Theta[n,n:1]*Theta[n,n:1]*v[1:n]) #nu(n)
	}
	return(list(theta=Theta, nu=v))
}	

ia = IA(N, gam)
theta.inf = ia$theta[N, 1:q];theta.inf
#[1] -0.1114082  0.2216191 -0.0309017 -0.3244026  0.3339867
nu.inf = ia$nu[N];nu.inf
#[1] 0.7816983

# c) Calculate and plot the autocovariance function defined by the parameters from b).
#    Compare with a).

gamma.ia <- sapply(h, FUN=gammaX, theta=theta.inf, sigma=sqrt(nu.inf));gamma.ia
sum(gamma.ia-gam)
#[1] 0.0261705 difference is less than 0.03
# smal difference in last term
plot(h, gamma.ia, type="h", col="blue", ylab = "ACVF") #from b)
lines(h+0.01, gam, type="h", col="red") #from a)

#d) If the two parameter sets are different then at most one of them 
# gives an invertible model. Which of them does that with probability one?

theta = runif(n=q, min = -2, max = 2);theta
Mod(polyroot(c(1,theta)))
#[1]  0.7795998  0.8878476  0.8878476  2.2131556 12.8982231 - non-invertible
Mod(polyroot(c(1,theta.inf)))
# [1] 1.438013 1.200501 1.438013 1.098225 1.098225 - invertible
#IA found the invertible model

# 5.2
# a)  For N = 1000, 
#     generate {Xt, t = 1, . . . , N} according to (5.1)MA(q)

q = 5
set.seed(1234)
theta = runif(q, -2, 2);theta
sigma2z = 1/(1+sum(theta^2));sigma2z

#(5.1) \ X_t = \sum_{j=0}^{q}\theta_j Z_{t-j}, \ \theta_0=1

N = 1000
sigmaz = sqrt(sigma2z)
#rlaplace -> lib(extraDistr)
library(extraDistr)
set.seed(1234)
Zt = rlaplace(n=N, mu=0, sigma=sigmaz/sqrt(2))
Xt = arima.sim(list(ma=theta), n=N, innov = Zt)
head(Xt)

#b) Calculate and plot the empirical autocovariance function
acf(Xt, type="covariance")

#c) Find {theta_inf,j | j = 1..q, v_inf} from the IA with input 
# {gamma(h) | h = 0..q} If the algortihm does not converge, 
# increase q to q' and use gamma(h) ;h = 0..q' as input
emp.gamma = as.vector(acf(Xt, type="covariance", plot=FALSE)$acf)
length(emp.gamma)
ia2 = IA(N, emp.gamma)
theta.inf2 = ia2$theta[N, 1:5];theta.inf2
nu.inf2 = ia2$nu[N];nu.inf2

# d) Compute and plot the autocovariance function defined the 
# parameters from the IA given in c). Compare with b)
gamma.acvf <- sapply(h, FUN=gammaX, theta=theta.inf2, sigma=sqrt(nu.inf2))
acf(Xt, type="covariance", lag.max = 15)
lines(h+0.08, gamma.acvf, type="h", col="green")
gamma.acvf
# [1]  1.00313296 -0.14888739  0.09920841  0.06194804 -0.24862657  0.23970684  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000
emp.gamma[1:(q+1)]
#[1]  1.00313296 -0.14888739  0.09920841  0.06194804 -0.24862657  0.22084301

# gamma.ia2 which is ACVF and  emp.gamma which is empirical ACVF 
# is the same expect last term. ACVF sligtly larger than empirical ACVF

# e) 
'''
The estimators theta_hat, sigma_hat from IA - satisfy the equation because, we see from the plot
that the emprical ACVF is approximately equal to the theoretical ACVF (mentioned in 5.1 a) 

What do you call these estimators?
Moment estimators.

How can you do this without using IA?
Could have used Newton-Rhapson/newton method - root finding algorithm
'''

# 5.4 An AR(6) model
# Do basic descriptive statistics and estimation for the there cases. 
# Present plots and for the estimation use both Yule Walker, least square and maximum likelihood
phi = c(0.40, 0.36, -0.47, 0.45, 0.21, -0.03)
set.seed(1234)

Mod(polyroot(c(1,phi)))

N1=100
Xt.1 = arima.sim(list(ar=phi), n=N1, innov = rnorm(N1) )
N2=1000
Xt.2 = arima.sim(list(ar=phi), n=N2, innov = rnorm(N2) )
N3=100000
Xt.3 = arima.sim(list(ar=phi), n=N3, innov = rnorm(N3) )

#descriptive statistics
summary(Xt.1);var(Xt.1)
summary(Xt.2);var(Xt.2)
summary(Xt.3);var(Xt.3)
par(mfrow=c(3,1))
plot(1:N1, Xt.1, type="l")
plot(1:N2, Xt.2, type="l")
plot(1:N3, Xt.3, type="l")
par(mfrow=c(3,1))
acf(Xt.1, type="covariance")
acf(Xt.2, type="covariance")
acf(Xt.3, type="covariance")


ar_order=6


ar.ols(x = Xt.1, aic = FALSE, order.max = 6)
ar(Xt.1, order=ar_order, method="yw")
ar(Xt.1, order=ar_order, method="mle")
ar(Xt.1, order=ar_order, method="burg")
'''
ols: 0.4659   0.2360  -0.4993   0.4478   0.1301   0.0372      
yw:	 0.4595   0.1501  -0.4686   0.4738  
mle: 0.5231   0.1386  -0.4991   0.5357  
burg:0.5001   0.1534  -0.5121   0.5110
'''

ar.ols(x = Xt.2, aic = FALSE, order.max = ar_order)
ar(Xt.2, order=ar_order, method="yw")
ar(Xt.2, order=ar_order, method="mle")
ar(Xt.2, order=ar_order, method="burg")
'''
ols: 0.3880   0.4106  -0.4998   0.4686   0.2161  -0.1100      
yw:	 0.3855   0.4031  -0.4921   0.4684   0.2103  -0.1042 
mle: 0.3855   0.4031  -0.4921   0.4684   0.2103  -0.1042 
burg:0.3868   0.4068  -0.4984   0.4713   0.2156  -0.1102 
'''

ar.ols(x = Xt.3, aic = FALSE, order.max = ar_order)
ar(Xt.3, order=ar_order, method="yw")
ar(Xt.3, order=ar_order, method="mle")
ar(Xt.3, order=ar_order, method="burg")
'''
ols: 0.3989   0.3646  -0.4747   0.4507   0.2132  -0.0366      
yw:	 0.3989   0.3646  -0.4746   0.4507   0.2131  -0.0366 
mle: 0.3988   0.3646  -0.4747   0.4506   0.2132  -0.0366 
burg:0.3988   0.3646  -0.4747   0.4506   0.2132  -0.0366 

phi = c(0.40, 0.36, -0.47, 0.45, 0.21, -0.03)
#Pretty close to ground truth when N is large
'''

spectrum(Xt.1, ..., )
par(mfrow = c(2,2))
spectrum(Xt.1 )  #method = c("pgram", "ar")
spectrum(Xt.1, spans = 3)
spectrum(Xt.1, spans = c(3,3))
spectrum(Xt.1, spans = c(3,5))

spectrum(ldeaths)
spectrum(ldeaths, spans = c(3,3))
spectrum(ldeaths, spans = c(3,5))
spectrum(ldeaths, spans = c(5,7))
spectrum(ldeaths, spans = c(5,7), log = "dB", ci = 0.8)

# acf, spectrum, spec.pgram, spectrum, arima, ar
# 1. Discuss which parameters that are significant for the different situations. 
# 2. How mabny observations are need for this size of the model? 
# 3. Comment a parametric versus a nonparametric estimate of the spectral density. 
# 4. Can you see any connection between the data and the empirical ACF or the spectral density?
# 5. You have used 3 different estimation methods for the autoregressive parameter vector. Do they differ much?

#Answere:
# 5: OLS, YW, MLE are equal for large N1
# 4: 
coeftest(arima(Xt.N1, order=c(p,0,0), method="ML"))
# 3:
# non-parametric: periodogram, parametric: AR estimation
# er to metoder i spectrum, en parametrisk og en ikke.
spectrum(Xt.N1, method="ar")
spectrum(Xt.N1, method="pgram")
#https://en.wikipedia.org/wiki/Spectral_density_estimation?fbclid=IwAR1azWt2zz6iM5O8gO5AaJpvHehBv28tM-dgPGbTGILagk5E5xX0zHy3kko #Techniques


