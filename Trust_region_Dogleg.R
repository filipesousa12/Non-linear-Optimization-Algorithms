library(numDeriv)
library(stats)


f<-function(x){
  exp(2*x[1])+x[2]*x[1]^2+x[2]^3+3*x[1]^2;
}

g<-function(x){
  grad(f,x)
}
h<-function(x){
  hessian(f,x)
}

xn<-c(1,2)
u<-1/4
eta<-3/4
dlt<-0.5
L<-10^(-4)

while(norm(g(xn),"2")>L){
  a<-t(g(xn))%*%h(xn)%*%g(xn)
  b<-norm(g(xn),"2")^2
  d<-dlt/norm(g(xn),"2")
  if (a>0){
    alpha<-min(d,b/a)
  } else {
    alpha<-d
  }
  q<-function(p){
    f(xn)+t(g(xn))%*%p+1/2*t(p)%*%h(xn)%*%p
  }
  cp<--alpha*g(xn)
  
  g1<-g(xn)
  h1<-h(xn)
  np<--solve(h1)%*%g1
  
  if(norm(cp,"2")==dlt){
    v<-cp
  } else if(norm(np,"2")<=dlt){
    v<-np
  } else {
    dogleg<-function(l){
    norm(np+l*(cp-np),"2")-dlt
  }
  l<-uniroot(dogleg,c(0,1))$root
  v<-l*cp+(1-l)*np
  }
  
  ro<-(f(xn)-f(xn+v))/(f(xn)-q(v))
  
  if (ro<=u){
    xv<-xn
  } else {
    xv<-xn
    xn<-xv+v
  }
  
  if(ro<=u){
    dlt<-dlt/2
  } else if(ro<=eta){
    dlt<-dlt
  } else {
    dlt<-2*dlt
  }
}
xn

