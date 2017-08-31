close all

n=80;
dx=2*pi/n;
x=0:dx:(2*pi-dx);
b=-cos(x)*dx*dx;
b=b';
sol=ConjugateGradient(b);
sol=sol-sum(sol)/n;
true=cos(x);
true=true';

figure,
plot(x,sol,'r',x,true,'g');
figure
plot(x,true-sol,'b')
max(true-sol)
