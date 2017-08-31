function Triangleshape(R)

str=['Triangle',num2str(R),'.run'];
Vf = fopen(str,'w');

H=1;
% R=L/H
%R=tan(pi/3);
L=R*H;

ns=258;
n=ns+1;

xs0=zeros(n,1);
ys0=zeros(n,1);
curv=zeros(n,1);

n1=80;
n2=(ns-n1)/2;
dh=H/n1;
dL=L/n2;
for i=1:n1/2+1
    xs0(i)=L/2;
    ys0(i)=(i-1)*dh;
end
for i=n1/2+2:ns/2+1;
    xs0(i)=L/2-(i-n1/2-1)*dL;
    ys0(i)=H/L/2*(xs0(i)+L/2);
end
for i=ns/2+2:2*n2+n1/2+1;
    xs0(i)=-L/2+(i-ns/2-1)*dL;
    ys0(i)=-H/L/2*(xs0(i)+L/2);
end
for i=2*n2+n1/2+1:n
    xs0(i)=L/2;
    ys0(i)=dh*(i-n);
end

figure,
plot(xs0,ys0);
axis([-L/2-0.5 L/2+0.5 -H-0.5 H+0.5]);
axis equal
fprintf(Vf,'xs0\t ys0\t curv\t \n');
for i=1:n
    fprintf(Vf,'%16.12f %16.12f %16.12f \n',xs0(i),ys0(i),curv(i));
end

fclose(Vf);

end
