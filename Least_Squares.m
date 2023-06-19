clear
t=[1,2,4,5,8];
y=[3.29,4.26,7.17,9.30,20.256];

function retval = fn(n , x)
  t=[1,2,4,5,8];
  y=[3.29,4.26,7.17,9.30,20.256];
  retval=(x(1)*exp(x(2)*t(n)))-y(n);
endfunction

function retval = F(x)
  F_1=[];
  for i=1:5
    F_1=[F_1;fn(i,x)];
   endfor
  retval=F_1;
endfunction

function retval = N_fn(n , x)
  t=[1,2,4,5,8];
  y=[3.29,4.26,7.17,9.30,20.256];
  retval=[exp(x(2)*t(n));x(1)*t(n)*exp(x(2)*t(n))];
endfunction

function retval = N_F(x)
  F_1=[];
  for i=1:5
    F_1=[F_1,N_fn(i,x)];
   endfor
  retval=F_1;
endfunction

N=10; #Numero iterações
L=0.0000001; #Norma erro
n=1;
u=10^(-4);
ro=1/2;
xn=[2.5;0.25];
while (n<=N & norm(N_F(xn)*F(xn))>L)
g=N_F(xn)*F(xn);
h=N_F(xn)*transpose(N_F(xn));
p=h\-g;
a=1;
%while(f(xn+a*p)>f(xn)+u*a*p.g1)
%a=*ro
%endwhile
xv=xn;
xn=xv+a*p;
n=n+1;
f=sum(F(xn));
endwhile
