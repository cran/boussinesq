

beq.lin.dimentionless <- function(t=0,x=seq(from=0,to=L,by=by),big=100000,by=L*0.01,L=1) {
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

beq.lin <-  function(t=0,x=seq(from=0,to=L,by=by),h1=1,h2=1,L=100,ks=0.01,q=ks,s=0.4,big=10^7,by=L/100,p=0.5) {
  
    D=(ks*(h1*p+h2*(1-p))/s) 
    
    out0 <- h2
    out1 <- (h1-h2)*beq.lin.dimentionless(t=(t*D)/L^2,x=x/L,big=100000)
    
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


