result
load xc.dat
load yc.dat

umax=-100;
ymax=-100;
vmax=-100;
xmax=-100;
vmin=100;
xmin=100;

n=129;
nmid=(n-1)/2+1;
for j=1:n
    if umax<uc(nmid,j)
        umax = uc(nmid,j);
        ymax = yc(j);
    end
end
umax
ymax

for i=1:n
    if vmax<vc(i,nmid)
        vmax = vc(i,nmid);
        xmax = xc(i);
    end
    if vmin>vc(i,nmid)
        vmin = vc(i,nmid);
        xmin = xc(i);
    end
end
vmax
xmax
vmin
xmin


