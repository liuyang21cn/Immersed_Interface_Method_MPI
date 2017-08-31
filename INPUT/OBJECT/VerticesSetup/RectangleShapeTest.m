clear all
close all
clc

a=0.75/2;
b=0.25/2;

ns=256;
n=256+1;
nb=10;db=b/nb;
na=ns/4-nb;da=(a+b)/na;

xs0=zeros(n,1);
ys0=zeros(n,1);
curv=zeros(n,1);

for i=1:nb
    xs0(i)=a+b;
    ys0(i)=(i-1)*db;
end

for i=nb+1:(2*na+nb)
    xs0(i)=(a+b)-(i-1-nb)*da;
    ys0(i)=b;
end

for i=(2*na+nb+1):(3*nb+2*na)
    xs0(i)=-(a+b);
    ys0(i)=b-(i-1-2*na-nb)*db;
end   
for i=(3*nb+2*na+1):(4*na+3*nb)
    xs0(i)=-(a+b)+(i-3*nb-2*na-1)*da;
    ys0(i)=-b;
end
for i=(4*na+3*nb+1):n
    xs0(i)=a+b;
    ys0(i)=-b+(i-4*na-3*nb-1)*db;
end

figure(1),
plot(xs0,ys0);
axis([-5 3 -5 3]);


t0=0;
tc=pi/0.4;
theta=0.75*pi+0.5*pi*(sin(0.8d0*0)*(1-exp(-t0/tc)));
ysc=1.25d0*(cos(0.8d0*t0)+1.0d0)*sin(pi/3.0d0);
xsc=1.25d0*(cos(0.8d0*t0)+1.0d0)*cos(pi/3.0d0);

thetat=0.25d0*pi*0.8d0*cos(0.8d0*t0)*(1.0d0-exp(-t0/tc))+ ...
    0.25d0*pi*(sin(0.8d0*t0))*(exp(-t0/tc)/tc);
xsct=1.25d0*(-0.8d0*sin(0.8d0*t0))*cos(pi/3.0d0);
ysct=1.25d0*(-0.8d0*sin(0.8d0*t0))*sin(pi/3.0d0);

xs1=zeros(n,1);
ys1=zeros(n,1);
for i=1:n
    xs1(i)=xsc+xs0(i)*cos(theta)-ys0(i)*sin(theta);
    ys1(i)=ysc+xs0(i)*sin(theta)+ys0(i)*cos(theta);
end
figure(2)
plot(xs1,ys1);
axis([-5 3 -5 3]);