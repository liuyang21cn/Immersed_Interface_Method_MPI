clear 
close all

img=imread('SMU_Mustang_Logo.svg.png');
bw=im2bw(img);
B=bwboundaries(bw);

A1=B(4);
A1=cell2mat(A1);
A1=A1*0.005;
x1=A1(:,1)*cos(-pi/2)-A1(:,2)*sin(-pi/2)-3;
y1=A1(:,1)*sin(-pi/2)+A1(:,2)*cos(-pi/2)+2;

A2=B(3);
A2=cell2mat(A2);
A2=A2*0.005;
x2=A2(:,1)*cos(-pi/2)-A2(:,2)*sin(-pi/2)-3;
y2=A2(:,1)*sin(-pi/2)+A2(:,2)*cos(-pi/2)+2;

A3=B(2);
A3=cell2mat(A3);
A3=A3*0.005;
x3=A3(:,1)*cos(-pi/2)-A3(:,2)*sin(-pi/2)-3;
y3=A3(:,1)*sin(-pi/2)+A3(:,2)*cos(-pi/2)+2;

figure,
hold on
plot(x1,y1,'.r');
plot(x2,y2,'.r');
plot(x3,y3,'.r');
hold off

n=length(x1);
curv=zeros(n,1);
Vf = fopen('Mustang.run','w');
fprintf(Vf,'xs0\t ys0\t curv\t \n');
for i=1:10:n
    fprintf(Vf,'%16.12f %16.12f %16.12f \n',x1(i),y1(i),curv(i));
end
fclose(Vf);
