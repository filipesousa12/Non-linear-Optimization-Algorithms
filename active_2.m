clear

A=[2,1;1,2];
b=[3;3];
R=[1,2];

xn = [5;0];


function retval = f(x)
  retval=0.5*x(1)^2+0.5*x(2)^2;
endfunction

function retval = g(x)
  retval=[x(1);x(2)];
endfunction


function retval = h(x)
  retval=[1,0;0,1];
endfunction


function retval=Rest_ativ(x)
  A=[2,1;1,2];
  b=[3;3];
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
  A=[2,1;1,2];
  b=[3;3];
  R=[1,2];
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
L = 10^-6;
S = 0;
X =[xn];
while (stop==0)
  Wi=Rest_ativ(xn);
  W=Wi;
  A_barra=A(Wi,:);
  n=rows(A_barra);
  if (columns(W)!=0)
    Z=null(A_barra);
  else
    Z=eye(2);
  endif
  A_r=A_barra\eye(n);
  if(columns(Z)==0)
    red_g=0;
  else
    red_g=transpose(Z)*g(xn);
  endif
  lambda=transpose(A_r)*g(xn);
  if(abs(red_g) < L)
    lambda=transpose(A_r)*g(xn);
     
    if(lambda>0)
      stop = 1
    else
      for (i=1:n)
          if(lambda(i,1)==min(lambda))
            W(i)=[];
            i=n;
          endif
      endfor
    endif
    
    A_barra=A(W,:);
    
    n=rows(A_barra);
    
    if (columns(W)!=0)
      Z=null(A_barra);
    else
      Z=eye(2);
    endif
    
   endif
   
    p=-Z*inverse(transpose(Z)*h(xn)*Z)*transpose(Z)*g(xn);
    alpha=min(1,a_otim(xn,p,W));
    if(alpha == 0)
      S = S +1
      if(S == 5)
        stop = 1
      endif
    endif
    xv=xn;
    xn=xv+alpha*p;
    X = [X , xn]
    N=N+1;
endwhile

