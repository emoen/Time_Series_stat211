library(itsmr)

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



# Numerically approximate causal AR(2) region
abline(a=1, b=0 , col=c("red"))
abline(a=-1, b=0 , col=c("red"))
abline(v=c(-2,2) , col=c("red"))

#Unifonr random sample in 2D
size <- 2000             #length of random number vectors
set.seed(1) 
x <- runif(size, min = -2, max = 2)          # generate samples from uniform distribution (0.0, 1.0)
y <-runif(size, min = -1, max = 1) 
df <-data.frame(x,y)

count=0
for( i in c(1:size) ) {
  a = specify(ar=c(x[i], y[i]))  
  phi = c(1,-a$phi)
  if (all(abs(polyroot(phi)) > 1)){
    points(x[i], y[i],  col = "red")
    count = count +1
  } else{
    points(x[i], y[i],  col = "gray")
  }
}
print(count/size)



#Unifonr random sample in 3D
size <- 200000             #length of random number vectors
set.seed(1) 
x <- runif(size, min = -2, max = 2)          # generate samples from uniform distribution (0.0, 1.0)
y <-runif(size, min = -1, max = 1)
z <-runif(size, min = -1, max = 1)
df <-data.frame(x,y)

xs = c()
ys = c()
zs = c()
count=0
for( i in c(1:size) ) {
  a = specify(ar=c(x[i], y[i], z[i]))  
  phi = c(1,-a$phi)
  if (all(abs(polyroot(phi)) > 1)){
    #points(x[i], y[i],  col = "red")
    xs = c(xs, x[i])
    ys = c(ys, y[i])
    zs = c(zs, z[i])
    count = count +1
  }
}
print(max(unlist(xs)))
bestx = which.max(unlist(xs))
print( paste(xs[bestx],  ys[bestx], zs[bestx], sep=" "))
print(max(unlist(ys)))
besty = which.max(unlist(ys))
print( paste(ys[bestz], ys[besty], zs[besty], sep=" "))
print(max(unlist(zs)))
bestz = which.max(unlist(zs))
print( paste(zs[bestz], ys[bestz], zs[bestz], sep=" "))
print(count)
print(count/size)


#Unifonr random sample in 4D
size <- 200000             #length of random number vectors
set.seed(1) 
x <- runif(size, min = -2, max = 2)          # generate samples from uniform distribution (0.0, 1.0)
y <-runif(size, min = -1, max = 1)
z <-runif(size, min = -1, max = 1)
v <-runif(size, min = -1, max = 1)

xs = c()
ys = c()
zs = c()
vs = c()
count=0
for( i in c(1:size) ) {
  a = specify(ar=c(x[i], y[i], z[i], v[i]))  
  phi = c(1,-a$phi)
  if (all(abs(polyroot(phi)) > 1)){
    #points(x[i], y[i],  col = "red")
    xs = c(xs, x[i])
    ys = c(ys, y[i])
    zs = c(zs, z[i])
    vs = c(vs, v[i])
    count = count +1
  }
}
print(max(unlist(xs)))
bestx = which.max(unlist(xs))
print( paste(xs[bestx],  ys[bestx], zs[bestx], vs[bestx], sep=" "))
print(max(unlist(ys)))
besty = which.max(unlist(ys))
print( paste(xs[bestz], ys[besty], zs[besty], vs[besty], sep=" "))
print(max(unlist(zs)))
bestz = which.max(unlist(zs))
print( paste(zs[bestz], ys[bestz], zs[bestz], vs[bestz], sep=" "))
print(max(unlist(vs)))
bestv = which.max(unlist(vs))
print( paste(xs[bestv], ys[bestv], zs[bestv], vs[bestv], sep=" "))
print(count)
print(count/size)

#Unifonr random sample in 5D
size <- 200000             #length of random number vectors
set.seed(1) 
x <- runif(size, min = -2, max = 2)          # generate samples from uniform distribution (0.0, 1.0)
y <-runif(size, min = -1, max = 1)
z <-runif(size, min = -1, max = 1)
v <-runif(size, min = -1, max = 1)
u <-runif(size, min = -1, max = 1)

xs = c()
ys = c()
zs = c()
vs = c()
us = c()
count=0
for( i in c(1:size) ) {
  a = specify(ar=c(x[i], y[i], z[i], v[i], u[i]))  
  phi = c(1,-a$phi)
  if (all(abs(polyroot(phi)) > 1)){
    #points(x[i], y[i],  col = "red")
    xs = c(xs, x[i])
    ys = c(ys, y[i])
    zs = c(zs, z[i])
    vs = c(vs, v[i])
    us = c(us, u[i])
    count = count +1
  }
}
print(max(unlist(xs)))
bestx = which.max(unlist(xs))
print( paste(xs[bestx],  ys[bestx], zs[bestx], vs[bestx], vs[bestx], us[bestx], sep=" "))
print(max(unlist(ys)))
besty = which.max(unlist(ys))
print( paste(xs[bestz], ys[besty], zs[besty], vs[besty], vs[besty], us[besty], sep=" "))
print(max(unlist(zs)))
bestz = which.max(unlist(zs))
print( paste(xs[bestz], ys[bestz], zs[bestz], vs[bestz], vs[bestz], us[bestz], sep=" "))
print(max(unlist(vs)))
bestv = which.max(unlist(vs))
print( paste(xs[bestv], ys[bestv], zs[bestv], vs[bestv], vs[bestv],  us[bestv], sep=" "))
print(max(unlist(us)))
bestu = which.max(unlist(us))
print( paste(xs[bestv], ys[bestv], zs[bestv], vs[bestv], vs[bestu], us[bestu], sep=" "))
print(count)
print(count/size)

#Unifonr random sample in 6D
size <- 200000             #length of random number vectors
set.seed(1) 
x <- runif(size, min = -2, max = 2)          # generate samples from uniform distribution (0.0, 1.0)
y <-runif(size, min = -1, max = 1)
z <-runif(size, min = -1, max = 1)
v <-runif(size, min = -1, max = 1)
u <-runif(size, min = -1, max = 1)
t <-runif(size, min = -1, max = 1)

xs = c()
ys = c()
zs = c()
vs = c()
us = c()
ts = c()
count=0
for( i in c(1:size) ) {
  a = specify(ar=c(x[i], y[i], z[i], v[i], u[i], t[i]))  
  phi = c(1,-a$phi)
  if (all(abs(polyroot(phi)) > 1)){
    #points(x[i], y[i],  col = "red")
    xs = c(xs, x[i])
    ys = c(ys, y[i])
    zs = c(zs, z[i])
    vs = c(vs, v[i])
    us = c(us, u[i])
    ts = c(ts, u[i])
    count = count +1
  }
}
print(max(unlist(xs)))
bestx = which.max(unlist(xs))
print( paste(xs[bestx],  ys[bestx], zs[bestx], vs[bestx], vs[bestx], us[bestx], sep=" "))
print(max(unlist(ys)))
besty = which.max(unlist(ys))
print( paste(xs[bestz], ys[besty], zs[besty], vs[besty], vs[besty], us[besty], sep=" "))
print(max(unlist(zs)))
bestz = which.max(unlist(zs))
print( paste(xs[bestz], ys[bestz], zs[bestz], vs[bestz], vs[bestz], us[bestz], sep=" "))
print(max(unlist(vs)))
bestv = which.max(unlist(vs))
print( paste(xs[bestv], ys[bestv], zs[bestv], vs[bestv], vs[bestv],  us[bestv], sep=" "))
print(max(unlist(us)))
bestu = which.max(unlist(us))
print( paste(xs[bestv], ys[bestv], zs[bestv], vs[bestv], vs[bestu], us[bestu], sep=" "))

print(max(unlist(ts)))
bestt = which.max(unlist(ts))
print( paste(xs[bestt], ys[bestt], zs[bestt], vs[bestt], vs[bestt], us[bestu], ts[bestt], sep=" "))
print(count)
print(count/size)

# Stupid f**king R method that doesnt return result but writes to output
checkLocal = function(a) {
  phi = c(1,-a$phi)
  theta = c(1,a$theta)
  if (all(abs(polyroot(phi)) > 1)){
    cat("Causal\n")
    return(1)
  }
  #else
  cat("Non-Causal\n")
  return(0)
  
  #if (all(abs(polyroot(theta)) > 1))
  #  cat("Invertible\n")
  #else
  #  cat("Non-Invertible\n")
}