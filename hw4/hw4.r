# 4.2 b) Use R and plot {ψj, j = 0, . . . , 50} when the parameteres are given by (4.12) (X_t= MA(INF))
#Consider a causal ARMA(2, 3)

library(itsmr)
phi = c(1.7,-.9)
theta = c(-1.4,.8,.1)
sigma2 = 1
#Using a function to check for invertibility and causality
check(list(phi=phi, theta=theta, sigma2 = sigma2))
## Causal
## Invertible
#Checking roots of AR-polynomial
Mod(polyroot(c(1,-phi)))
## [1] 1.054093 1.054093
#Checking roots of MA-polynomial
Mod(polyroot(c(1,theta)))
## [1] 1.022124 1.022124 9.571780


psi = numeric(51)
psi[1] = 1
psi[2] = phi[1]+theta[1]
psi[3] = phi[1]*psi[2]+phi[2]*psi[1]+theta[2]
psi[4] = phi[1]*psi[3]+phi[2]*psi[2]+theta[3]
for(k in 5:51){
    psi[k] = phi[1]*psi[k-1]+phi[2]*psi[k-2]
}
plot(0:50,psi, xlab = "k", ylab = expression(psi[k]), type= "b")


#4.3 c) 
#Implement the results in R, compute and plotf
#(h); h = 0; : : : ; 50g with parameter
#values given by (3). Check your calculations with help of an R-function.

m = 50
M = matrix(0, ncol=m+1,nrow=m+1)
diag(M) = 1
M[2,2] = 1-phi[2]; M[1,2:3] = -phi; M[m+1,m] = -phi[1]
for(i in 1:(m-1)){
    M[seq(i+1,i+2),i] = -phi
}

M[1:5,1:5]
## [,1] [,2] [,3] [,4] [,5]
## [1,] 1.0 -1.7 0.9 0.0 0
## [2,] -1.7 1.9 0.0 0.0 0
## [3,] 0.9 -1.7 1.0 0.0 0
## [4,] 0.0 0.9 -1.7 1.0 0
## [5,] 0.0 0.0 0.9 -1.7 1
M[(m-3):(m+1),(m-3):(m+1)]
## [,1] [,2] [,3] [,4] [,5]
## [1,] 1.0 0.0 0.0 0.0 0
## [2,] -1.7 1.0 0.0 0.0 0
## [3,] 0.9 -1.7 1.0 0.0 0
## [4,] 0.0 0.9 -1.7 1.0 0
## [5,] 0.0 0.0 0.9 -1.7 1
s = c(1+theta[1]*psi[2]+theta[2]*psi[3]+theta[3]*psi[4],
theta[1]+theta[2]*psi[2]+theta[3]*psi[3],
theta[2]+theta[3]*psi[2],theta[3],rep(0,nrow(M)-4))
gamma.h = solve(M,s)
plot(0:m,gamma.h, col = 4, type = "h", xlab= "h", ylab = "Autocovariance")

#4.5
# Use R and draw the rectangle defined by (4.11) in a φ1φ2-coordinate system
phi1 = seq(from = -2.5, to = 2.5, length = 51) 
plot(phi1,1+phi1,lty="dashed",type="l",xlab="",ylab="",cex.axis=.8,ylim=c(-1.5,1.5))
abline(a = -1, b = 0, lty="dashed")
abline(a = 1, b = -1, lty="dashed")
title(ylab=expression(phi[2]),xlab=expression(phi[1]),cex.lab=.8)
polygon(x = phi1[6:46], y = 1-abs(phi1[6:46]), col="gray")
lines(phi1,-phi1^2/4)
text(0,-.5,expression(phi[2]<phi[1]^2/4),cex=.7)
text(1.2,.5,expression(phi[2]>1-phi[1]),cex=.7)
text(-1.75,.5,expression(phi[2]>1+phi[1]),cex=.7)