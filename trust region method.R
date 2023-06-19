library(numDeriv)
f <- function(x){
  exp(2*x[1])+x[2]*x[1]^2+x[2]^3+3*x[1]^2
}
f(xk)
xk <- matrix(c(1,2),ncol=1,nrow = 2)
d_k <- 0.5
mu<- 1/4
eta<-3/4

for(i in 1:100){
  g = grad(f,xk)
  h = hessian(f,xk)
  a = t(g)%*%h%*%g
  b = norm(g,"2")^2
  d = d_k/norm(g,"2")
  
  alfa = d
  
  if(a>0){
    alfa = min(d,b/a)
  }

  
  pk = -alfa*g
  
  qk = f(xk)+t(g) %*% pk+0.5*t(pk) %*% h %*% pk
  
  ro <- (f(xk)-f(xk+pk))/(f(xk)-qk)
  
  if(ro > mu){
    xk = xk+pk
  }
  else{
    xk=xk
  }
  
  if(ro<=mu){
    d_k=0.5*d_k
  }
  
  else if(mu<= ro && ro<eta){
    d_k = d_k
    
  }
  
  else{
    d_k = 2*d_k
  }
  
  if(d_k >3){
    d_k=3
  }
  
  
  
  
}
xk
f(xk)
