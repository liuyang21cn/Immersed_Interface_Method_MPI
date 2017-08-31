imax=0;
jmax=0;
mm=-1000;

data=rhsp;
data0=rhsp0;

for i=1:240
    for j=1:240
        if(abs(data(i,j)-data0(i,j))>mm)
            mm=abs(data(i,j)-data0(i,j));
            imax=i;
            jmax=j;
        end
    end
end
imax
jmax
mm