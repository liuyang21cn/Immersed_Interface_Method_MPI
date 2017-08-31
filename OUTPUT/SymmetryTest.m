m=2401;
n=961;
result;


% rhsp
rhsp1 = rhsp(1:m,1:(n-1)/2);
rhsp2 = rhsp(1:m,1:(n-1)/2);

for i = 1:m
    for j = 1:(n-1)/2
        rhsp2(i,j)=rhsp(i,(n+1)-j);
    end
end
figure,
surf(rhsp1-rhsp2)
title('rhsp1-rhsp2')
fprintf('rhsp \t')
location(m,(n-1)/2,abs(rhsp1),abs(rhsp2));

% p
p=p';
p1 = p(1:m,1:(n-1)/2);
p2 = p(1:m,1:(n-1)/2);
for i = 1:m
    for j = 1:(n-1)/2
        p2(i,j)=p(i,(n+1)-j);
    end
end
figure,
surf(p1-p2),
title('p')
fprintf('p \t')
location(m,(n-1)/2,p1,p2);

% uc
uc=uc';
uc1 = uc(1:m,1:(n-1)/2);
uc2 = uc(1:m,1:(n-1)/2);

for i = 1:m
    for j = 1:(n-1)/2
        uc2(i,j)=uc(i,(n+1)-j);
    end
end
figure,
surf(uc1-uc2),
title('uc')
fprintf('uc \t')
location(m,(n-1)/2,uc1,uc2);


% vc
vc=vc';
vc1 = vc(1:m,1:(n-1)/2);
vc2 = vc(1:m,1:(n-1)/2);
for i = 1:m
    for j = 1:(n-1)/2
        vc2(i,j)=vc(i,(n+1)-j);
    end
end
figure,
% surf(vc1-vc2)
surf(abs(vc1)-abs(vc2)),
title('vc')
fprintf('vc \t')
location(m,(n-1)/2,abs(abs(vc1)),abs(abs(vc2)));

