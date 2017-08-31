imax=0;
jmax=0;
mm=-1000;

data=d;

for i=1:240
    for j=1:240
        if(max(d(i,j))>mm)
            mm=d(i,j);
            imax=i;
            jmax=j;
        end
    end
end
imax
jmax
mm