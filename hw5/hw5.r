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
#Library for rlaplace: extraDistr
#mu, sigma location and scale parameters. Scale must be positive.
library(extraDistr)
set.seed(1234)
Zt = rlaplace(n=N, mu=0, sigma=sigmaz/sqrt(2))
Xt = arima.sim(list(ma=theta), n=N, innov = Zt)
head(Xt)

#b) Calculate and plot the empirical autocovariance function
acf(Xt, type="covariance")

#c) Find {theta_inf,j | j = 1..q, v_inf} from the IA with input 
# {gamm(h) | h = 0..q} If the algortihm does not converge, 
# increase q to q' and use gamma(h) ;h = 0..q' as input
gamma.emp = as.vector(acf(Xt, type="covariance", plot=FALSE)$acf)
length(gamma.emp)
ia2 = IA(N, gamma.emp)
theta.inf2 = ia2$theta[N, 1:5];theta.inf2
nu.inf2 = ia2$nu[N];nu.inf2
