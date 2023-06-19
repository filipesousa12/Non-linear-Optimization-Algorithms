  library(numDeriv)
f <- function(x){
  x[2]*exp(x[1])+x[1]^3-3*x[2] 
}

xk <- matrix(c(1,1),ncol=1,nrow = 2)

B<- diag(nrow(xk))
f(xk)
for(i in 1:2){
  g = grad(f,xk)
  g
  p = -solve(B)%*%g
  a = 1

  xv <- xk
  xk <- xk+a*p
  xk
  sk <- xk-xv
  yk <- grad(f,xk)-grad(f,xv)
  f(xk)
  
  B <- B-((B%*%sk)%*%t(B%*%sk))/as.numeric(t(sk)%*%B%*%sk)+(yk%*%t(yk))/as.numeric(t(yk)%*%sk)
  B
  
  
}

round(xk,4)
round(f(xk),4)

norm(g,"2")
