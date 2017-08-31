sp=zeros(241,241);
sp0=zeros(241,241);

for i=1:241
    for j=1:241
        sp(i,j)=ux(i,j)*vy(i,j);
        sp0(i,j)=ux0(i,j)*vy0(i,j);
    end
end

surf(sp-sp0)