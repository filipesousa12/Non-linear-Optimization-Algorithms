clear
function retval=f(x,n)
  if(n==0)
    retval=x(1)^2+x(2)^2;
  elseif (n==1)
    retval=[2*x(1);2*x(2)];
  elseif (n==2)
    retval=[2,0;0,2];
   endif
endfunction

function retval=g(x,n)
  if(n==0)
    retval=x(1)-1;
  elseif (n==1)
    retval=[1;0];
  elseif (n==2)
    retval=[0,0;0,0];
   endif
endfunction

xn=[1.001; 0.001];
N=0;
mu=10^(-4);

while(N<1)
lambda=mu/g(xn,0);
A=g(xn,1);
H=f(xn,2)-lambda*g(xn,2);

Z=null(A');
n=rows(A);
Ar=(A\eye(n))';

nab_B=f(xn,1)-mu/g(xn,0)*g(xn,1);

D=lambda^2;
p1=-Z*inverse(Z'*H*Z)*Z'*nab_B;
lam=Ar'*(H*p1+nab_B);
p2=(Z*inverse(Z'*H*Z)*Z'*H-eye(n))*Ar*inverse(D)*lam;

p=p1+mu*p2;
alpha=1;
while (g(xn+alpha*p,0)<0)
    alpha=alpha/2;
endwhile

xv=xn;
xn=xv+alpha*p;

mu=mu/1.1;
N=N+1;
endwhile
