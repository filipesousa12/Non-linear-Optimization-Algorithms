library(numDeriv)
k<-2
f<-function(x){
  4*(x[1]-2)^2+(x[2]-1)^2
}

f2<-function(x){
  grad(f,x)
}
f3<-function(x){
  hessian(f,x)
}

g<-function(x){
  1-x[1]^2+-x[2]^2
}

g2<-function(x){
  grad(g,x)
}

g3<-function(x){
  hessian(g,x)
}

Lag<-function(x,l){
  f(x)-t(l)%*%g(x)
}

Lag2<-function(x,l){
  f2(x)-t(l)%*%g2(x)
}

Lag3<-function(x,l){
  f3(x)-l*g3(x)
}

xn<-c(2,1)
ln<-0


N<-1000 #Numero iterações
L<-10^(-4)
n<-0

while (n<=(N-1) & norm(g(xn),"2")>L){
  h1<-matrix(0,k+1,k+1)
  h1[1:2,1:2]<-Lag3(xn,ln)
  h1[1:2,3]<--g2(xn)
  h1[3,1:2]<--g2(xn)
  g1<-c(Lag2(xn,ln),-g(xn))
  g1
  p<--solve(h1)%*%g1
  p
  p1<-p[1:2]
  p1
  v<-p[3]
  xv<-xn;
  xn<-xv+p1;
  xn
  lv<-ln
  ln<-ln+v
  ln
  n<-n+1;
}
h1
