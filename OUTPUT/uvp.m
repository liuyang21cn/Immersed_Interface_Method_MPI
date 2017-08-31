function uvp(m,n)

%clear
close all

fid=fopen('p.dat');
p=fread(fid,[m n],'double');
p=p';

fid=fopen('uc.dat');
uc=fread(fid,[m n],'double');
uc=uc';

fid=fopen('vc.dat');
vc=fread(fid,[m n],'double');
vc=vc';
fclose('all');

figure(1),
surf(p)
figure(2),
surf(uc)
figure(3)
surf(vc)
%figure(4),
%surf(rhsp)

end