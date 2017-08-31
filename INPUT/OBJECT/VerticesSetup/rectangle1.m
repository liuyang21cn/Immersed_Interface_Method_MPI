

str = sprintf('vobj%i',index);
Vf = fopen(str,'w');

a=0.75/2;
b=0.25/2;

ns=256;
n=256+1;
nb=10;db=b/nb;
na=ns/4-nb;da=(a+b)/na;

xs0=zeros(n,1);
ys0=zeros(n,1);
curv=zeros(n,1);
alfa=zeros(n,1);

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

plot(xs0,ys0);
axis([-1 1 -1 1]);
fprintf(Vf,'xs0\t ys0\t curv\t \n');
for i=1:n
    fprintf(Vf,'%16.12f %16.12f %16.12f \n',xs0(i),ys0(i),curv(i));
end

fclose(Vf);
