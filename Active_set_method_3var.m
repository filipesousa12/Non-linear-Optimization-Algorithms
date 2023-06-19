clear
A=[2,1,1;1,2,1;1,1,2];
b=[2;3;4];
R=[1,2,3];

xn=[0;1;2];

X=xn

function retval=f(x)
  retval=1/2*x(1)^2+1/2*x(2)^2+1/2*x(3)^2;
endfunction

function retval=g(x)
  retval=[x(1);x(2);x(3)];
endfunction

function retval=h(x)
  retval=[1,0,0;0,1,0;0,0,1];
endfunction

function retval=Rest_ativ(x)
  A=[2,1,1;1,2,1;1,1,2];
b=[2;3;4];
  n=rows(A);
  W=[];
  j=0;
  for (i=1:n)
    if (A(i,:)*x==b(i))
      j=j+1;
      W(j)=i;
    endif
  endfor
  retval=W;
endfunction

function retval=a_otim(x,p,W)
  A=[2,1,1;1,2,1;1,1,2];
b=[2;3;4];
  R=[1,2,3];
  minimo=10^6;
  S=setdiff(R,W);
  m=length(S);
  for (i=1:m)
    if (-(A(S(i),:)*x-b(S(i)))/(A(S(i),:)*p)<minimo && A(S(i),:)*p<0)
      minimo=-(A(S(i),:)*x-b(S(i)))/(A(S(i),:)*p);
    endif
  endfor
  retval=minimo
endfunction

stop=0;
N=0;
L=10^(-6);
while (stop==0 && N<=100)
  Wi=Rest_ativ(xn);
  W=Wi;
  A_barra=A(Wi,:);
  n=rows(A_barra);
  if (columns(W)!=0)
    Z=null(A_barra);
  else
    Z=eye(3);
  endif
  A_r=A_barra\eye(n);
  if(columns(Z)==0)
    red_g=0;
  else
    red_g=transpose(Z)*g(xn);
  endif
  while((abs(red_g)<L||columns(Z)==0) && stop==0)
  
    if(columns(Wi)==0)
      stop=1
    else
      lambda=transpose(A_r)*g(xn);
    endif
    
    if (min(lambda)<0)
        for (i=1:n)
          if(lambda(i,1)==min(lambda))
            W(i)=[];
            i=n;
          endif
        endfor
     else
        stop=1
    endif
    A_barra=A(W,:);
    n=rows(A_barra);
    if (columns(W)!=0)
      Z=null(A_barra);
    else
      Z=eye(2);
    endif
    A_r=A_barra\eye(n);
    red_g=transpose(Z)*g(xn);
  endwhile
  
  p=-Z*inverse(transpose(Z)*h(xn)*Z)*transpose(Z)*g(xn);
  alpha=min(1,a_otim(xn,p,W));
  xv=xn;
  xn=xv+alpha*p;
  N=N+1;
  X=[X,xn];
 endwhile
 xn
 f(xn)