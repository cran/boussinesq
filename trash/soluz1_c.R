
#
# Author: Emanuele Cordano
#
# Example: analytic solution to boussinesq equation with two-dounded Dirichlet condition
#

steadystatesol <- function(x=seq(from=0,to=L,by=by),h1=1,h2=1,L=100,ks=0.01,q=ks,by=L*0.01) {

  beta0=h1^2
  beta1=(h2^2-h1^2)/L+q/ks*L 
  
  out=(-q/ks*x^2+beta1*x+beta0)^0.5
  
  return(out)
}

wsol_a <- function(t=0,x=seq(from=0,to=L,by=by),big=100000,by=L*0.01,L=1) {
#
# t is dimensionless and resacaled with L^2/D  and q*L^2/
# x is dimensionless and resacaled with L
# coefficient is dimnsionless and rescaled with (q*L^2)/(D*s)
#

 sum=0
 for (n in 0:big) sum=sum+1/((2*n+1)*pi)^3*(1-exp(-(2*n+1)^2*pi^2*t))*sin((2*n+1)*x*pi)

  return(sum)

}

approxsol <-  function(t=0,x=seq(from=0,to=L,by=by),h1=1,h2=1,L=100,ks=0.01,q=ks,s=0.4,big=10^7,by=L/100) {
  
    D=(ks*(h1+h2))/(2*s) 
    
    out0 <- steadystatesol(x=x,h1=h1,h2=h2,L=L,ks=ks,q=0)
    out1 <- q*L^2/(D*s)*wsol_a(t=(t*D)/L^2,x=x/L,big=100000)
    
    return(out0+out1)
}


wsol_a_p2 <- function(t=0,x=seq(from=0,to=L,by=by),big=100000,by=L*0.01,L=1) {
#
# t is dimensionless and resacaled with L^2/D  and q*L^2/
# x is dimensionless and resacaled with L
# coefficient is dimnsionless and rescaled with (q*L^2)/(D*s)
#

 sum=0
 for (n in 1:big) sum=sum-2/(pi*n)*exp(-n^2*pi^2*t)*sin(n*pi*x) 
 #sum=sum+1/((2*n+1)*pi)^3*(1-exp(-(2*n+1)^2*pi^2*t))*sin((2*n+1)*x*pi)

  sum=sum+1-x
  
  return(sum)

}

approxsol_p2 <-  function(t=0,x=seq(from=0,to=L,by=by),h1=1,h2=1,L=100,ks=0.01,q=ks,s=0.4,big=10^7,by=L/100,p=0.5) {
  
    D=(ks*(h1*p+h2*(1-p))/s) 
    
    out0 <- h2
    out1 <- (h1-h2)*wsol_a_p2(t=(t*D)/L^2,x=x/L,big=100000)
    
    return(out0+out1)
}





# 
# plot_approxsol <- function(t=0,x=seq(from=0,to=L,by=by),h1=1,h2=1,L=100,ks=0.01,q=ks,s=0.4,big=10^7,by=L/100,xlab="x [m]",ylab=" h [m]",main="bla") {
#  
#   x <- seq(from=0,to=L,by=by)
#   
#   h <- array(NA,c(length(t),length(x)))
#   
#   for (i in 1:length(t)) {
#      h[i,] <- approxsol(t=t[i],x=x,h1=h1,h2=h2,L=L,ks=ks,q=q,s=s,big=big,by=by)
#   }
#   
#  #to do  plot(x,h[i,],xlab=xlab,ylab=ylab,main=main,type="l")
#   
#  
# }
# 

plot_approxsol_p2 <- function(t=10^(0:6),x=seq(from=0,to=L,by=by),h1=2,h2=1,L=100,ks=0.01,q=ks,s=0.4,big=10^7,by=L/100,xlab="x [m]",ylab=" h [m]",main="Water table depth vs x",loc="topright") {
 
  x <- seq(from=0,to=L,by=by)
  
  h <- array(NA,c(length(t),length(x)))
  
  for (i in 1:length(t)) {
  
     h[i,] <- approxsol_p2(t=t[i],x=x,h1=h1,h2=h2,L=L,ks=ks,q=q,s=s,big=big,by=by)
  }
  
  plot(x,h[1,],xlab=xlab,ylab=ylab,main=main,type="l",ylim=c(min(h),max(h)),lty=1)
  
   for (i in 2:length(t)) {
      lines(x,h[i,],lty=i)

   }
  
  legend(loc,lty=1:length(t),legend=t)
  
}

# asymtotic solution

asymptotic_solution <-  function(x=seq(from=0,to=L,by=by),h1=1,h2=1,L=100,ks=0.01,q=0,s=0.4,big=10^7,by=L/100,p=0.5) {
  
#    D=(ks*(h1*p+h2*(1-p))/s) 
    b0 <- h1^2
    b1 <- (h2^2-h1^2)/L
    
    out <- (-q/ks*x^2+b1*x+b0)^0.5
    
  #  out0 <- h2
  #  out1 <- (h1-h2)*wsol_a_p2(t=(t*D)/L^2,x=x/L,big=100000)
    
    return(out)
}


