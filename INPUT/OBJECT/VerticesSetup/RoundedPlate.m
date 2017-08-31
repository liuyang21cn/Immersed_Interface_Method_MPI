clear all
close all
clc

Vf = fopen('RoundedPlate.run','w');

a=0.75;
b=0.25;

ns=256;
n=256+1;

xs0=zeros(n,1);
ys0=zeros(n,1);
curv=zeros(n,1);

dalfa=2*pi/ns;

xf=a/2;
yf=b/2;
ssl=(2*a+pi*b)/ns;
angle=-ssl/yf;
xsi=-xf;
ysi=yf;

for m=0:ns
    if ((xsi+ssl)<=xf)&&(xsi>=-xf)&&(ysi>0)
        xsi=xsi+ssl;
    elseif ((xsi-ssl)>=-xf)&&(xsi<=xf)&&(ysi<0)
        xsi=xsi-ssl;
    elseif (xsi>=xf)&&(abs(atan2(ysi,(xsi-xf))+angle)<=pi/2)
        tmp=(xsi-xf)*cos(angle)-ysi*sin(angle);
        ysi=(xsi-xf)*sin(angle)+ysi*cos(angle);
        xsi=tmp+xf;
    elseif (xsi<=-xf)&&(abs(atan2(ysi,(xsi+xf))+angle)>=pi/2)
        tmp=(xsi+xf)*cos(angle)-ysi*sin(angle);
        ysi=(xsi+xf)*sin(angle)+ysi*cos(angle);
        xsi=tmp-xf;
    elseif ((xsi+ssl)>xf)&&(xsi<xf)&&(ysi>0)
        tangle=-(xsi+ssl-xf)/yf;
        xsi=xf;
        tmp=(xsi-xf)*cos(tangle)-ysi*sin(tangle);
        ysi=(xsi-xf)*sin(tangle)+ysi*cos(tangle);
        xsi=tmp+xf;
    elseif ((xsi-ssl)<-xf)&&(xsi>-xf)&&(ysi<0)
        tangle=(xsi-ssl+xf)/yf;
        xsi=-xf;
        tmp=(xsi+xf)*cos(tangle)-ysi*sin(tangle);
        ysi=(xsi+xf)*sin(tangle)+ysi*cos(tangle);
        xsi=tmp-xf;
    elseif (xsi>=xf) 
        tsl=yf*(abs(atan2(ysi,(xsi-xf))+angle)-pi/2);
        ysi=-yf;
        xsi=xf-tsl;
    elseif (xsi<=-xf)
        tsl=yf*(pi/2-abs(atan2(ysi,(xsi+xf))+angle));
        ysi=yf;
        xsi=-xf+tsl;
    end
    xs0(n-m)=xsi;
    ys0(n-m)=ysi;
end

temp1 = sqrt(a*a+b*b)/2;
for m=1:n
    temp2=sqrt(xs0(m)*xs0(m)+ys0(m)*ys0(m));
    if(temp2<temp1) 
        curv(m)=0;
    else
        curv(m)=2/b;
    end
end

figure,
plot(xs0,ys0);
axis([-1 1 -1 1]);
fprintf(Vf,'xs0\t ys0\t curv\t \n');
for i=1:n
    fprintf(Vf,'%16.12f %16.12f %16.12f \n',xs0(i),ys0(i),curv(i));
end

fclose(Vf);