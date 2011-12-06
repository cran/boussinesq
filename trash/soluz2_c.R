
#
# Author: Emanuele Cordano
#
# Example: Song'and Lockington's analytic solution to boussinesq equation with two-dounded Dirichlet condition
#

song_so <- function (t=0.5,x=1.0,s=0.4,h1=1,ks=0.01,nmax=4,alpha=1) {
 
    xi <- x*(2*s*(alpha+1)/(h1*ks*t^(alpha+1)))^0.5
    ax <- build_a(n=nmax,lambda=alpha/(alpha+1))
    xi0 <- sum(ax)^(-0.5)
    out <- h1*t^alpha*song_hval(xi=xi,xi0=xi0,a=ax*xi0^2)
 
    return(out)
}

song_hval <- function(xi,xi0,a) {
  

  temp <- 1-xi/xi0
  temp[temp<0] <- 0
  
  a0 <- 0
  out <- array(a0,length(temp))
  
  for (n in 1:length(a)) out <- out+a[n]*temp^n
  return(out)
}

build_a <- function(n=4,lambda=0) {
  a <- array(NA,n)
  a[1]=1/4
  a[2]=(2*lambda-1)/16
  for (i in 3:n) {
    a[i]=(2*lambda+1-i)/i^2*a[i-1]
    for (k in 2:(i-1)) a[i]=a[i]-2*(i+1)/i*a[k]*a[i+1-k]
 
    
  } 
  return(a)
}  

