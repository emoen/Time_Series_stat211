#or AR(3) the area of causal roots it is less than 30%, for AR(4) it is less than 17%, for AR(5) it is less than 9% and for #AR(6) it is less than 5%

size <- 200000           #Uniform random sample in 6D  
x <- runif(size, min = -2, max = 2) 
y <-runif(size, min = -1, max = 1)
z <-runif(size, min = -1, max = 1)
v <-runif(size, min = -1, max = 1)
u <-runif(size, min = -1, max = 1)
t <-runif(size, min = -1, max = 1)

count=0
for( i in c(1:size) ) {
  a = specify(ar=c(x[i], y[i], z[i], v[i], u[i], t[i]))  
  phi = c(1,-a$phi)
  if (all(abs(polyroot(phi)) > 1)){
    count = count +1
  }
}
print(count/size)





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
