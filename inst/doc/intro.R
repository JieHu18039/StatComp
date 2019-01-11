## ------------------------------------------------------------------------
set.seed(1)
x<-c(0:4) ; p<-c(0.1,0.2,0.2,0.2,0.3)
cp<-cumsum(p); n<-1000
r<-x[findInterval(runif(n),cp)+1] #has already get zhe sample of size 1000 by incerse transform
set.seed(1)
y<-sample(x,size=n,replace = TRUE,prob=p)
counts <- table(c(rep(0,n),rep(1,n)),c(r, y))
barplot(counts, main="The distribution of X",
        xlab="X",ylab='Count', col=c("darkblue","red"),
        legend = c('Inverse','sample'), beside=TRUE) #compare result

## ------------------------------------------------------------------------
r.beta<-function(a,b,size){
  k<-0; j<-0; y<-numeric(size)
  c<-dbeta((a-1)/(a+b-2),a,b) #The max value of beta(a,b)'s df
while (k < size) {
  u <- runif(1)
  j <- j + 1
  x <- runif(1) #random variate from g
  row.x<-dbeta(x,a,b)/c
  if ( row.x> u) {
    #we accept x
    k <- k + 1
    y[k] <- x
  }
}
  acp<-k/j #The accept rate
  return(list(acp,y))
}

x<-r.beta(3,2,1000)
x[[1]] #The accept rate
hist(x[[2]],breaks=seq(0,1,0.05),freq = FALSE,main = "Histogram of beta(3,2)",xlab = "beta(3,2)")
y<-seq(0,1,0.01)
lines(y,dbeta(y,3,2),col="red")

## ------------------------------------------------------------------------
n<-1000;r<-4;beta<-2
lamada<-rgamma(n,r,beta)
y<-rexp(n,lamada)
hist(y,breaks = seq(0,max(y)+1,0.1),freq = F)
f1<-function(r,beta,x){
  r*(beta^r)/((beta+x)^(r+1))
} #The pdf of Pareto distribution 
x<-seq(0,max(y)+1,0.01)
lines(x,f1(r,beta,x),col="red")

## ------------------------------------------------------------------------
set.seed(12)
beta_33<-function(x,m){
  n<-length(x)
  theta<-numeric(n)
  for(i in 1:n){
  y<-runif(m,0,x[i]) #sample Y
  theta[i]<-mean(x[i]*(y^2)*((1-y)^2)/beta(3,3))
}
  return(theta)
}

x<-seq(0.1,0.9,0.1)
m<-100000

x1<-beta_33(x,m)
x1<-round(x1,5) #estime result
x2<-pbeta(x,3,3)
a<-rbind(x,x1,x2)
rownames(a)<-c("X","MC","pbeta")
colnames(a)<-c("x1","x2","x3","x4","x5","x6","x7","x8","x9")
knitr::kable(a,format="markdown",caption = "Comparation of them",align = "c")

## ------------------------------------------------------------------------
 RL<-function(sigma,m,antithetic = TRUE){
  y<-runif(m)
  if(antithetic) v<-1-y else v<-runif(m)
  x1<-sigma*sqrt(-2*log(y))
  x2<-sigma*sqrt(-2*log(v))
  return(cbind(x1,x2))
 }
 sigma<-2 #the value of sigma
 m<-1000
 mc1<-RL(sigma,m) # antithetic samples
 mc2<-RL(sigma,m,antithetic = FALSE) #independent samples
 sd1<-sd(rowMeans(mc1))  
 sd2<-sd(rowMeans(mc2))  
 print(c(sd1,sd2,sd1/sd2))

## ------------------------------------------------------------------------
 g<-function(x) (x^2)*(exp((-x^2)/2))/sqrt(2*3.14159)
 f1<-function(x) exp(-(x-1))
 f2<-function(x) 1/x^2
 x<-seq(1,10,.01)
 gs<-c(expression(g(x)==x^2*e^{-x^2/2}/sqrt(2*pi)),expression(f1(x)==e^-(x-1)),expression(f2(x)==1/x^2))
 
 plot(x,g(x), type = "l", ylab = "",ylim = c(0,1), lwd = 0.25,col=1,main='(A)')
 lines(x, f1(x), lty = 2, lwd = 0.25,col=2)
 lines(x, f2(x), lty = 3, lwd = 0.25,col=3)
 legend("topright", legend = gs, lty = 1:3, lwd = 0.25, inset = 0.01,col=1:3)

## ------------------------------------------------------------------------
 gs1<-c(expression(g^2/f1),expression(g^2/f2))
 plot(x,g(x)^2/f1(x),type="l", ylab = "",ylim = c(0,0.3),lwd = 0.25,col=1,main='(B)')
 lines(x,g(x)^2/f2(x), lty = 2, lwd = 0.25,col=2)
 legend("topright", legend = gs1,lty = 1:2, lwd = 0.25, inset=0.1,col=1:2) 

## ------------------------------------------------------------------------

 m<-1000
 n<-10000
 result1<-numeric(m)
 for(i in 1:m){
 y1<-rexp(n,1)+1
 result1[i]<-mean(g(y1)/f1(y1))
 }
 m1<-mean(result1)
 s1<-sd(result1)
 
 result2<-numeric(m)
 for(i in 1:m){
   u<-runif(n)
   y2<-1/(1-u)
   result2[i]<-mean(g(y2)/f2(y2))
 } 
 m2<-mean(result2)
 s2<-sd(result2)
 b<-matrix(c(m1,s1,m2,s2),2,2)
 colnames(b)<-c("f1","f2")
 rownames(b)<-c("theta","se")
 knitr::kable(b,format="markdown",caption = "Comparation of them",align = "c")

## ------------------------------------------------------------------------
G<-function(x){
  n<-length(x)
  a<-seq(1-n,n-1,2)
  x.i<-sort(x)
  G.hat<-sum(a*x.i)/(n*n*mean(x))
  return(G.hat)
} # you can estimate a G.hat if there comes a sample

set.seed(1012)
# if X is standard lognormal
n=500
m<-1000
G.hat1<-numeric(m)
for(i in 1:m){
  y<-rnorm(n)
  x<-exp(y) # then x is standard lognormal
  G.hat1[i]<-G(x)
}
result1<-c(mean(G.hat1),quantile(G.hat1,probs=c(0.5,0.1)))
names(result1)<-c("mean","median","deciles")
print(result1)
hist(G.hat1,breaks=seq(min(G.hat1)-0.01,max(G.hat1)+0.01,0.01),freq =F,main = "Histogram of G",xlab = "standard lognormal")

# if X is uniform
G.hat2<-numeric(m)
for(i in 1:m){
  x<-runif(n) # then x is uniform
  G.hat2[i]<-G(x)
}
result2<-c(mean(G.hat2),quantile(G.hat2,probs=c(0.5,0.1)))
names(result2)<-c("mean","median","deciles")
print(result2)
hist(G.hat2,breaks =seq(min(G.hat2)-0.01,max(G.hat2)+0.01,0.01) ,freq =F,main = "Histogram of G",xlab = "uniform")

#if x is  Bernoulli(0.1)
G.hat3<-numeric(m)
for(i in 1:m){
  x<-rbinom(n,1,0.1) # then x is Bernoulli(0.1)
  G.hat3[i]<-G(x)
}
result3<-c(mean(G.hat3),quantile(G.hat3,probs=c(0.5,0.1)))
names(result3)<-c("mean","median","deciles")
print(result3)
hist(G.hat3,breaks=seq(min(G.hat3)-0.01,max(G.hat3)+0.01,0.01),freq =F,main = "Histogram of G",xlab = "Bernoulli(0.1)")


## ------------------------------------------------------------------------
n<-500
m<-1000
mu<-1
sigma<-1
G1<-function(n,mu,sigma){
  y<-rnorm(n,mu,sigma)
  x<-exp(y)
  G.sample<-G(x)
  return(G.sample)
}
G.sp<-numeric(m)
for(i in 1:m){
G.sp[i]<-G1(n,mu,sigma)
}
#the  approximate 95% conffidence interval is
CI<-c(mean(G.sp)-sd(G.sp)*qt(0.975,m-1),mean(G.sp)+sd(G.sp)*qt(0.975,m-1))
print(CI)
cover.rate<-sum(I(G.sp>CI[1]&G.sp<CI[2]))/m
print(cover.rate)

## ------------------------------------------------------------------------
library(MASS)
m<-1000
sigma1<-matrix(c(1,0.2,0.2,1),2,2)
p.value1<-p.value2<-p.value3<-numeric(m)
for(i in 1:1000){
x1<-mvrnorm(20,mu=c(0,0),Sigma = sigma1)
p.value1[i]<-cor.test(x1[,1],x1[,2],method = "pearson")$p.value
p.value2[i]<-cor.test(x1[,1],x1[,2],method = "spearman")$p.value
p.value3[i]<-cor.test(x1[,1],x1[,2],method = "kendall")$p.value
}
#the power
power<-c(mean(p.value1<=0.05),mean(p.value2<=0.05),mean(p.value3<=0.05))
names(power)<-c("pearson","spearman","kendall")
print(power)

## ----echo=FALSE----------------------------------------------------------
library(boot)
library(bootstrap)
library(DAAG)

## ------------------------------------------------------------------------
set.seed(1)
LSAT<-c(576, 635, 558, 578,666, 580, 555, 661, 651, 605, 653, 575, 545, 572, 594)
GPA<-c(339, 330, 281, 303, 344, 307, 300, 343, 336, 313, 312, 274, 276, 288, 296)
x<-cbind(LSAT,GPA)
n <-length(GPA)
R.hat <- cor(x[,1],x[,2]) 
R.jack <- numeric(n) 
for(i in 1:n)
  { R.jack[i] <-cor(x[-i,1],x[-i,2]) } 
bias.jack <- (n-1)*(mean(R.jack)-R.hat) 
se.jack <- sqrt((n-1)*mean((R.jack-R.hat)^2)) 
round(c(original=R.hat,bias=bias.jack, se=se.jack),3)

## ------------------------------------------------------------------------
y<-c(3,5,7,18,43,85,91,98,100,130,230,487)
boot.mean <- function(y,i) mean(y[i]) 
de <- boot(data=y,statistic=boot.mean, R =1000) 
ci <- boot.ci(de,type=c("norm","basic","perc","bca")) 
ci$norm[2:3]
ci$basic[4:5]
ci$percent[4:5]
ci$bca[4:5]

## ------------------------------------------------------------------------
attach(scor)
sigma.hat<-cov(scor)
n<-nrow(scor)
value<-eigen(sigma.hat)$value #get all Eigenvalues
theta.hat<-value[1]/sum(value)
theta.jack<- numeric(n) 
for(i in 1:n){ 
  sigma.1<-cov(scor[-i,])
  value.1<-eigen(sigma.1)$value
  theta.jack[i] <-value.1[1]/sum(value.1) 
  } 
bias.jack <- (n-1)*(mean(theta.jack)-theta.hat) 
se.jack <- sqrt((n-1)*mean((theta.jack-theta.hat)^2)) 
round(c(original=theta.hat,bias=bias.jack, se=se.jack),3)

## ------------------------------------------------------------------------
attach(ironslag)
n <- length(magnetic) #in DAAG ironslag

e1 <- e2 <- e3 <- e4<- matrix(0,n*(n-1)/2,2)

# for n*(n-1)/2-fold cross validation
# fit models on leave-two-out samples 
i=0 
for (k in 1:(n-1)) { 
  for (j in (k+1):n) {
  i=i+1  
  
  y <- magnetic[-c(k,j)] 
  x <- chemical[-c(k,j)]
  
 J1 <- lm(y ~ x) 
 yhat1 <- J1$coef[1] + J1$coef[2] * chemical[c(k,j)] 
 e1[i,] <- magnetic[c(k,j)] - yhat1
 
 J2 <- lm(y ~ x + I(x^2)) 
 yhat2 <- J2$coef[1] + J2$coef[2] * chemical[c(k,j)] + J2$coef[3] * chemical[c(k,j)]^2 
 e2[i,] <- magnetic[c(k,j)] - yhat2
 
 J3 <- lm(log(y) ~ x) 
 logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[c(k,j)] 
 yhat3 <- exp(logyhat3) 
 e3[i,] <- magnetic[c(k,j)] - yhat3
 
 J4 <- lm(log(y) ~ log(x)) 
 logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[c(k,j)]) 
 yhat4 <- exp(logyhat4) 
 e4[i,] <- magnetic[c(k,j)] - yhat4
 }
}

c(mean(e1^2),mean(e2^2),mean(e3^2),mean(e4^2))

## ----echo=FALSE----------------------------------------------------------
library(survival)
library(mvtnorm)
library(gam)
library(splines)
library(foreach)
library(RANN)
library(energy)
library(boot)
library(Ball)


## ------------------------------------------------------------------------
set.seed(1)
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"])) 
y <- sort(as.vector(weight[feed == "linseed"])) 
detach(chickwts)

CramerVonMisesTwoSamples <- function(x1,x2){
  Fx1<-ecdf(x1)
  Fx2<-ecdf(x2)
  n<-length(x1)
  m<-length(x2)
  w1<-sum((Fx1(x1)-Fx2(x1))^2)+sum((Fx1(x2)-Fx2(x2))^2) 
  w2<-w1*m*n/((m+n)^2)
  return(w2)
}  #get the Cramer-von Mises statistic

R <- 999
z <- c(x, y)
K<- 1:26
n<-length(x)
reps <- numeric(R)
t0 <-  CramerVonMisesTwoSamples (x,y)
for (i in 1:R) {
  k<- sample(K, size = n, replace = FALSE) 
  x1 <- z[k]
  y1 <- z[-k] #complement of x1 
  reps[i] <- CramerVonMisesTwoSamples (x1,y1)
  } 
  p <- mean(abs(c(t0, reps)) >= abs(t0)) 
  round(p,3)
  
  hist(reps, main = "", freq = FALSE, xlab = "Cramer-von Mises statistic", breaks = "scott")

  points(t0, 0, cex = 1, pch = 16) #observed T


## ------------------------------------------------------------------------
m<-1e2;k<-3;p<-2;mu<-0.5;R<-999;
p.values <- matrix(0,m,3)

Tn<-function(z,ix,sizes,k){
n1<-sizes[1];n2<-sizes[2];n<-n1+n2
if(is.vector(z)) z<-data.frame(z,0)
z<-z[ix, ]
NN<-nn2(data=z,k=k+1)
block1<-NN$nn.idx[1:n1,-1] 
block2<-NN$nn.idx[(n1+1):n,-1] 
i1<-sum(block1 <= n1)
i2<-sum(block2 >n1) 
(i1+i2)/(k*n)
}

eqdist.nn<-function(z,sizes,k){
boot.obj<-boot(data=z,statistic=Tn,R=R,sim="permutation",sizes=sizes,k=k)
ts<-c(boot.obj$t0,boot.obj$t)
p.value<-mean(ts>=ts[1])
list(statistic=ts[1],p.value=p.value)
}

## ------------------------------------------------------------------------
n1<-n2<-50;n<-n1+n2;N=c(n1,n2)  
for(i in 1:m){
x<-rnorm(n1,0,1)
y<-rnorm(n2,0,2)
z<-c(x,y)
p.values[i,1]<-eqdist.nn(z,N,k)$p.value
p.values[i,2]<-eqdist.etest(z,sizes=N,R=R)$p.value
p.values[i,3]<-bd.test(x=x,y=y,R=999,seed=i*12345)$p.value
}
alpha<-0.1
pow<-colMeans(p.values<alpha)
pow

## ------------------------------------------------------------------------
n1<-n2<-50;n<-n1+n2;N=c(n1,n2)  
for(i in 1:m){
x<-rnorm(n1,0,1)
y<-rnorm(n2,1,2)
z<-c(x,y)
p.values[i,1]<-eqdist.nn(z,N,k)$p.value
p.values[i,2]<-eqdist.etest(z,sizes=N,R=R)$p.value
p.values[i,3]<-bd.test(x=x,y=y,R=999,seed=i*12345)$p.value
}
alpha<-0.1
pow<-colMeans(p.values<alpha)
pow

## ------------------------------------------------------------------------
n1<-n2<-50;n<-n1+n2;N=c(n1,n2)  
for(i in 1:m){
x<-rt(n1,1)
y<-c(rnorm(n2/2,5,1),rnorm(n2/2))
z<-c(x,y)
p.values[i,1]<-eqdist.nn(z,N,k)$p.value
p.values[i,2]<-eqdist.etest(z,sizes=N,R=R)$p.value
p.values[i,3]<-bd.test(x=x,y=y,R=999,seed=i*12345)$p.value
}
alpha<-0.1
pow<-colMeans(p.values<alpha)
pow

## ------------------------------------------------------------------------
n1<-10;n2<-100;n<-n1+n2;N=c(n1,n2) 
for(i in 1:m){
x<-rnorm(10)
y<-rnorm(100)
z<-c(x,y)
p.values[i,1]<-eqdist.nn(z,N,k)$p.value
p.values[i,2]<-eqdist.etest(z,sizes=N,R=R)$p.value
p.values[i,3]<-bd.test(x=x,y=y,R=999,seed=i*12345)$p.value
}
alpha<-0.1
pow<-colMeans(p.values<alpha) 
pow

## ------------------------------------------------------------------------
 f1<-function(x,u,lamada){
    pi<-3.14159265
    return(lamada/(pi*(lamada^2+(x-u)^2)))
} 

cauchy.chain<-function(n,sigma,x0,u,lamada){
    x<-NULL
    x[1]<-x0
    e=runif(n)
    k<-0
    for(i in 2:n){
      y<-rnorm(1,x[i-1],sigma)
      if(e[i]<=(f1(y,u,lamada)/f1(x[i-1],u,lamada)))
      {x[i]<-y}
      else
      {x[i]<-x[i-1]
      k<-k+1}
    }
    return(list(x=x,k=k))
  } 
  n<-10000
  cauchy<-cauchy.chain(n,sigma=0.5,x0=0,u=0,lamada=1)
  refuse.pr<-cauchy$k/n
  refuse.pr
  #qq plot
  par(mfrow=c(1,1))
  qqplot(rcauchy(n-1000,0,1),cauchy$x[1001:n])
  qqline(cauchy$x[1001:n])
  #hist
  hist(cauchy$x[1001:n],freq=F,main="cauchy",breaks=60)
  curve(f1(x,0,1),add=TRUE)

## ------------------------------------------------------------------------
  w <-0.25 #width of the uniform support set 
  m <- 5000 #length of the chain 
  burn <- 1000 #burn-in time 
  y<-c(125,18,20,34)
  x <- numeric(m) #the chain
  
  prob <- function(b, y) { 
    # computes (without the constant) the target density 
    if (b < 0 || b >= 1) return (0)
    return((1/2+b/4)^y[1] * ((1-b)/4)^y[2] * ((1-b)/4)^y[3] * (b/4)^y[4])
  }
  
  u <- runif(m) #for accept/reject step
  v <- runif(m, -w, w) #proposal distribution 
  x[1] <-0.25 
  for (i in 2:m) { 
    z <- x[i-1] + v[i] 
    if (u[i] <= prob(z,y) / prob(x[i-1],y)) 
      x[i] <-z 
    else 
      x[i] <- x[i-1] 
  }
  
   xb <- x[(burn+1):m]
   xc<-mean(xb)
   print(xc)  #estimation value of theta
   print(y/sum(y))
   print(c(1/2+xc/4,(1-xc)/4,(1-xc)/4,xc/4))

## ------------------------------------------------------------------------
set.seed(1)
S.k<-function(a,k)  pt(a^2*(k-1)/(k-a^2),df=k-1,lower.tail=F)-pt(a^2*k/(k+1-a^2),df=k,lower.tail=F)
k<-c(4:25,100,500,1000)
root<-numeric(length(k))
j<-0
for(i in k){
j<-j+1  
solution<-uniroot(S.k,c(1e-5,sqrt(i)-(1e-5)),k=i)
root[j]<-solution$root
}
s<-cbind(k,root)
knitr::kable(s,format="markdown",caption = "Comparation of them",align = "c")

## ------------------------------------------------------------------------
set.seed(1)
Gelman.Rubin <- function(psi) { 
  # psi[i,j] is the statistic psi(X[i,1:j]) 
  # for chain in i-th row of X 
  psi <- as.matrix(psi) 
  n <- ncol(psi) 
  k <- nrow(psi)
  psi.means <- rowMeans(psi) #row means
  B <- n * var(psi.means) #between variance est. 
  psi.w <- apply(psi, 1, "var") #within variances 
  W <- mean(psi.w) #within est. 
  v.hat <- W*(n-1)/n + (B/n) #upper variance est. 
  r.hat <- v.hat / W #G-R statistic 
  return(r.hat) 
  }


prob <- function(b, y) { 
  # computes (without the constant) the target density 
  if (b < 0 || b >= 1) return (0)
  return((1/2+b/4)^y[1] * ((1-b)/4)^y[2] * ((1-b)/4)^y[3] * (b/4)^y[4])
}

theta.chain<-function(m,y,x1){
w <-0.25 #width of the uniform support set 
u <- runif(m) #for accept/reject step
v <- runif(m, -w, w) #proposal distribution 

x <- numeric(m) #the chain
x[1] <-x1 
for (i in 2:m) { 
  z <- x[i-1] + v[i] 
  if (u[i] <= prob(z,y) / prob(x[i-1],y)) 
    x[i] <-z 
  else 
    x[i] <- x[i-1] 
}
return(x)
}


b<- 1000 #burn-in time 
k<-4     #number of chains to generate
m <- 15000 #length of the chain 
y<-c(125,18,20,34)
x1 <- c(0.1,0.35,0.6,0.85)

#generate the chains 
x<- matrix(0, nrow=k, ncol=m) 
for (i in 1:k) 
  x[i, ] <- theta.chain(m,y,x1[i])

#compute diagnostic statistics 
psi <- t(apply(x, 1, cumsum)) 
for (i in 1:nrow(psi)) 
  psi[i,] <- psi[i,] / (1:ncol(psi)) 
print(Gelman.Rubin(psi))

#plot psi for the four chains 
par(mfrow=c(2,2),mar=c(2,2,2,2)) 
for (i in 1:k) 
  plot(psi[i, (b+1):m], type="l", xlab=i, ylab=bquote(psi)) 
par(mfrow=c(1,1)) #restore default

#plot the sequence of R-hat statistics 
rhat <- rep(0, m) 
for (j in (b+1):m) 
  rhat[j] <- Gelman.Rubin(psi[,1:j]) 
plot(rhat[(b+1):m], type="l", xlab="", ylab="R") 
abline(h=1.1, lty=2)

## ----echo=FALSE----------------------------------------------------------
library(nloptr)

## ------------------------------------------------------------------------
f<-function(y,sita,eta){
#sita mean scale parameter
#eta mean the location parameter
1/(sita*3.141592653*(1+((y-eta)/sita)^2))
}
# the cauchy pdf

pdf<-function(x,sita,eta,lower.tail=TRUE){
 if(lower.tail) res<-integrate(f,lower = -Inf,upper = x,rel.tol=.Machine$double.eps^0.25,sita=sita,eta=eta)
 else res<-integrate(f,lower = x,upper = Inf,rel.tol=.Machine$double.eps^0.25,sita=sita,eta=eta)
  return(res$value)
}
pdf(x=0,sita = 1,eta = 0)
pcauchy(0,location = 0,scale = 1)

pdf(x=2,sita = 2,eta =1,lower.tail = F )
pcauchy(2,location = 1,scale = 2,lower.tail = F)

## ----echo=FALSE----------------------------------------------------------
        dat <- rbind(Genotype=c('AA','BB','OO','AO','BO','AB'),
                     Frequency=c('p2','q2','r2','2pr','2qr','2pq',1),
                     Count=c('nAA','nBB','nOO','nAO','nBO','nAB','n'))
    knitr::kable(dat,format='markdown',caption = "Comparation of them",align = "c")

## ------------------------------------------------------------------------

# Mle function
eval_f0 <- function(x,x1,n.A=28,n.B=24,nOO=41,nAB=70) {
  #x[1] mean p , x1[1] mean p0
  #x[2] mean q , x1[2] mean q0
  r1<-1-sum(x1)
  nAA<-n.A*x1[1]^2/(x1[1]^2+2*x1[1]*r1)
  nBB<-n.B*x1[2]^2/(x1[2]^2+2*x1[2]*r1)
  r<-1-sum(x)
  return(-2*nAA*log(x[1])-2*nBB*log(x[2])-2*nOO*log(r)-
           (n.A-nAA)*log(2*x[1]*r)-(n.B-nBB)*log(2*x[2]*r)-nAB*log(2*x[1]*x[2]))
}


# constraint function 
eval_g0 <- function(x,x1,n.A=28,n.B=24,nOO=41,nAB=70) {
  return(sum(x)-0.999999)
}

opts <- list("algorithm"="NLOPT_LN_COBYLA",
             "xtol_rel"=1.0e-8)
mle<-NULL
r<-matrix(0,1,2)
r<-rbind(r,c(0.2,0.35))# the beginning value of p0 and q0
j<-2
while (sum(abs(r[j,]-r[j-1,]))>1e-8) {
res <- nloptr( x0=c(0.3,0.25),
               eval_f=eval_f0,
               lb = c(0,0), ub = c(1,1), 
               eval_g_ineq = eval_g0, 
               opts = opts, x1=r[j,],n.A=28,n.B=24,nOO=41,nAB=70 )
j<-j+1
r<-rbind(r,res$solution)
mle<-c(mle,eval_f0(x=r[j,],x1=r[j-1,]))
}
r  #the result of EM algorithm
mle #the max likelihood values


## ------------------------------------------------------------------------
set.seed(1)
attach(mtcars)
formulas <- list( mpg ~ disp, mpg ~ I(1 / disp), mpg ~ disp + wt, mpg ~ I(1 / disp) + wt )
#for loops
out <- vector("list", length(formulas))
for (i in seq_along(formulas)) { out[[i]] <-lm(formulas[[i]]) }
out

#lapply
lapply(formulas,lm)

## ------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ] })

#for loops
for(i in seq_along(bootstraps)){
  print(lm(mpg~disp,data =bootstraps[[i]]))
}

#lapply
lapply(bootstraps,lm,formula=mpg~disp)


## ------------------------------------------------------------------------
#in exercise 3
rsq <- function(mod) summary.lm(mod)$r.squared
#for loops
for (i in seq_along(formulas)) {
 print( rsq(lm(formulas[[i]])))
  }
#lapply
lapply(lapply(formulas,lm),rsq)
#in exercise 4
#for loops
for(i in seq_along(bootstraps)){
  print(rsq(lm(mpg~disp,data =bootstraps[[i]])))
}

#lapply
lapply(lapply(bootstraps,lm,formula=mpg~disp),rsq)

## ------------------------------------------------------------------------
#using anonymous function
trials <- replicate( 100, t.test(rpois(10, 10), rpois(7, 10)), simplify = FALSE )
p_value<-function(mod) mod$p.value
sapply(trials, p_value)
#using [[ directly


## ------------------------------------------------------------------------
x <- list(a = c(1:3), b = c(4:8))

lapply.f<-function(x,f,...){
   r<-Map(f,x,...)
   n<-length(r[[1]])
   return(vapply(r,as.vector,numeric(n)))
}
lapply.f(x,mean)
lapply.f(x,quantile)

## ----echo=F--------------------------------------------------------------
library(microbenchmark)

## ------------------------------------------------------------------------
set.seed(1)
my.chisq.statistic<-function(x,y){
  if(!is.vector(x) && !is.vector(y))
    stop("at least one of 'x' and 'y' is not a vector")
  if(typeof(x)=="character" || typeof(y)=="character")
    stop("at least one of 'x' and 'y' is not a numeric vector")
  if(any(x<0) || anyNA(x)) 
    stop("all entries of 'x' must be nonnegative and finite")
  if(any(y<0) || anyNA(y)) 
    stop("all entries of 'y' must be nonnegative and finite")
  if((n<-sum(x))==0) 
    stop("at least one entry of 'x' must be positive")
  if((n<-sum(x))==0) 
    stop("at least one entry of 'x' must be positive")
  if(length(x)!=length(y)) 
    stop("'x' and 'y' must have the same length")
  DNAME<-paste(deparse(substitute(x)),"and",deparse(substitute(y)))
  METHOD<-"Pearson's Chi-squared test"
  x<-rbind(x,y)
  nr<-as.integer(nrow(x));nc<-as.integer(ncol(x))
  sr<-rowSums(x);sc<-colSums(x);n<-sum(x)
  E<-outer(sr,sc,"*")/n
  STATISTIC<-sum((x - E)^2/E)
  names(STATISTIC)<-"X-squared"
  structure(list(statistic=STATISTIC,method=METHOD,data.name=DNAME),class="htest")
}

x<-sample(100:200,1e4,replace = TRUE)
y<-sample(100:200,1e4,replace = TRUE)
m<-rbind(x,y)
my.chisq.statistic(x,y)$statistic
chisq.test(m)$statistic
ts<-microbenchmark(mean1=my.chisq.statistic(x,y)$statistic,meanR2=chisq.test(m)$statistic)
summary(ts)[,c(1,3,5,6)]

## ------------------------------------------------------------------------
my.table<-function(...,dnn = list.names(...),deparse.level = 1){
    list.names <- function(...) {
        l <- as.list(substitute(list(...)))[-1L]
        nm <- names(l)
        fixup <- if (is.null(nm)) 
            seq_along(l)
        else nm == ""
        dep <- vapply(l[fixup], function(x) switch(deparse.level + 
            1, "", if (is.symbol(x)) as.character(x) else "", 
            deparse(x, nlines = 1)[1L]), "")
        if (is.null(nm)) 
            dep
        else {
            nm[fixup] <- dep
            nm
        }
    }
    args <- list(...)
    if (!length(args)) 
        stop("nothing to tabulate")
    if (length(args) == 1L && is.list(args[[1L]])) {
        args <- args[[1L]]
        if (length(dnn) != length(args)) 
            dnn <- if (!is.null(argn <- names(args))) 
                argn
            else paste(dnn[1L], seq_along(args), sep = ".")
    }
    bin <- 0L
    lens <- NULL
    dims <- integer()
    pd <- 1L
    dn <- NULL
    for (a in args) {
        if (is.null(lens)) 
            lens <- length(a)
        else if (length(a) != lens) 
            stop("all arguments must have the same length")
        fact.a <- is.factor(a)
        if (!fact.a) {
            a0 <- a
            a <- factor(a)
        }
        ll <- levels(a)
        a <- as.integer(a)
        nl <- length(ll)
        dims <- c(dims, nl)
        dn <- c(dn, list(ll))
        bin <- bin + pd * (a - 1L)
        pd <- pd * nl
    }
    names(dn) <- dnn
    bin <- bin[!is.na(bin)]
    if (length(bin)) 
        bin <- bin + 1L
    y <- array(tabulate(bin, pd), dims, dimnames = dn)
    class(y) <- "table"
    y
}

mya<-myb<-c(1,seq(1,4))            #example
my.table(mya,myb)
table(mya,myb)
microbenchmark(t1=my.table(mya,myb),t2=table(mya,myb)) 

