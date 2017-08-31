% clear
close all

m=1201;
n=451;

fid=fopen('iop.dat');
iop=fread(fid,[m n],'int16');

fid=fopen('iou.dat');
iou=fread(fid,[m n],'int16');

fid=fopen('iov.dat');
iov=fread(fid,[m n],'int16');
fclose('all');

loadObj
load xc.dat
load yc.dat
load xf.dat
load yf.dat

% iop
figure,
hold on
% plot(obj0(:,1),obj0(:,2),'k-')
plot(obj1(:,1),obj1(:,2),'k-')
% plot(obj2(:,1),obj2(:,2),'k-')
% plot(obj3(:,1),obj3(:,2),'k-')
% plot(obj4(:,1),obj4(:,2),'k-')
for i=1:m
    for j=1:n
        if(iop(i,j)==257)
            plot(xc(i),yc(j),'.r');
        elseif(iop(i,j)==258)
            plot(xc(i),yc(j),'.k');
        elseif(iop(i,j)==180)
            plot(xc(i),yc(j),'.b');
        elseif(iop(i,j)==140)
            plot(xc(i),yc(j),'.g');
        end
    end
end
hold off
title('iop');
xlabel('x');ylabel('y');
% axis([-5 3 -5 3])
axis equal

% iou
figure,
hold on
plot(obj1(:,1),obj1(:,2),'k-')
% plot(obj2(:,1),obj2(:,2),'k-')
% plot(obj3(:,1),obj3(:,2),'k-')
% plot(obj4(:,1),obj4(:,2),'k-')
for i=1:m
    for j=1:n
        if(iou(i,j)==257)
            plot(xf(i+1),yc(j),'.r');
        elseif(iou(i,j)==258)
            plot(xf(i+1),yc(j),'.k');
        elseif(iou(i,j)==180)
            plot(xf(i+1),yc(j),'.b');
        elseif(iou(i,j)==140)
            plot(xf(i+1),yc(j),'.g');
        end
    end
end
hold off
% axis([-2 2 -2 2])
title('iou');
xlabel('x');ylabel('y');
axis equal

% iov
figure,
hold on
plot(obj1(:,1),obj1(:,2),'k-')
% plot(obj2(:,1),obj2(:,2),'k-')
% plot(obj3(:,1),obj3(:,2),'k-')
% plot(obj4(:,1),obj4(:,2),'k-')
for i=1:m
    for j=1:n
        if(iov(i,j)==257)
            plot(xc(i),yf(j+1),'.r');
        elseif(iov(i,j)==258)
            plot(xc(i),yf(j+1),'.k');
        elseif(iov(i,j)==180)
            plot(xc(i),yf(j+1),'.b');
        elseif(iov(i,j)==140)
            plot(xc(i),yf(j+1),'.g');
        end
    end
end
hold off
% axis([-2 2 -2 2])
title('iov');
xlabel('x');ylabel('y');
axis equal
