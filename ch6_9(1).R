####   例 6.5.5 独立抽样--混合正态模型
rm(list=ls())
m <- 10000  #length of chain
xt <- numeric(m)
p <- .2                #mixing parameter
n <- 30               #sample size
mu <- c(0, 5)      #parameters of the normal densities
sigma <- c(1, 1)

sigma_p = c(0.5,1,2,4,.8) #variance of the proposal distribution

# generate the observed sample
i <- sample(1:2, size=n, replace=TRUE, prob=c(p, 1-p))
x <- rnorm(n, mu[i], sigma[i])

# generate the independence sampler chain
u <- runif(m)
xt[1] <- .5
k=0
for (i in 2:m) {
  y <- rnorm(1, xt[i-1], sigma_p[3])
  fy <- y * dnorm(x, mu[1], sigma[1]) +
    (1-y) * dnorm(x, mu[2], sigma[2])
  fx <- xt[i-1] * dnorm(x, mu[1], sigma[1]) +
    (1-xt[i-1]) * dnorm(x, mu[2], sigma[2])
  
  r <-  exp((y-xt[i-1])^2 / 2*(sigma_p[3])) * prod(fy / fx)
  
  if (u[i] <= r) xt[i] <- y else
    xt[i] <- xt[i-1]
  k=k+1
}

# plot for convergence diagnostic purpose
par(mfrow=c(1,2))
plot(xt, type="l", ylab="p")
hist(xt[1001:m], main="", xlab="p", prob=TRUE)

z = seq(min(xt[1001:m]),max(xt[1001:m]),length=1000)
lines(z,dnorm(z,mean(xt[1001:m]),sd(xt[1001:m])))

print(mean(xt[1001:m]))
print(var(xt[1001:m]))
print(k)/m