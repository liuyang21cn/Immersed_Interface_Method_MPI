close all

load fbypp.dat
[m,n]=size(fbypp);

figure,
% for i=1:m
%     plot(fbypp(i,1),fbypp(i,2),'r');
% end
hold on
plot(fbypp(:,1),fbypp(:,2),'r')
plot(fbypp(:,1),fbypp(:,3),'b')
hold off


