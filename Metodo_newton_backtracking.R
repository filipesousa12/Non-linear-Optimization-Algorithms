library(numDeriv)
f <- function(x){
  (x[1]-1)^2+x[2]^2  
}

xk <- matrix(c(10,-0),ncol=1,nrow = 2)

for(i in 1:100){
  g = grad(f,xk)
  h = hessian(f,xk)
  p = -solve(h)%*%g
  a = 1
  while (f(xk+a*p)> f(xk) - 0.5*a*t(p)%*%grad(f,xk)){
    
    a = a*(1/2)
  }
  xk <- xk+a*p
  
}

round(xk,4)
round(f(xk),4)

norm(g,"2")
    