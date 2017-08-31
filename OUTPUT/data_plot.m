clear all

result

load xc.dat
load yc.dat

ms=1024;

kip=1;

icontour=1;
iplot=1;
ishape=1;
iquiver=0;
iplotall=0;

igiven=0;
label=0;
level=1000;

mc=size(xc,1);
nc=size(yc,1);

ixcs=1;
ixce=mc;
jycs=1;
jyce=nc;

wov=[-2.5:0.1:2.5];
phv1=[-8:0.2:8];
phv2=[-0.02:0.002:0.02];
phv=[phv1,phv2];
pv=[-4:0.2:4];

figure(1)
if icontour==1
    if iplot==1
        if igiven==0
            H=pcolor(xc,yc,wo);
            shading interp;
            caxis([-5 5]);
            axis equal;
        else
            cs=contour(xc,yc,wo,wov);
        end
    end
    if iplot==2
        if igiven==0
            cs=contour(xc,yc,ph,level);
        else
            cs=contour(xc,yc,ph,phv);
        end
    end
    if iplot==3
        if igiven==0
            H=pcolor(xc,yc,p);
            shading interp;
            caxis([-1 1]);
            axis equal;
            %      [cs,h]=contour(xc,yc,p,level);
            %      surf(xc,yc,p)
            %      set(h,'ShowText','on','TextStep',get(h,'LevelStep')*4)
        else
            cs=contour(xc,yc,p,pv);
        end
    end
    if iplot==4
        %    H=pcolor(xc,yc,d);
        %    shading interp;
        %    caxis([-1 1]);
        %    axis equal;
        cs=contour(xc,yc,d,level);
    end
    if label==1
        clabel(cs)
    end
    hold on
end
if iquiver==1
    quiver(xc(ixcs:kip:ixce),yc(jycs:kip:jyce),...
        uc(jycs:kip:jyce,ixcs:kip:ixce),vc(jycs:kip:jyce,ixcs:kip:ixce))
    hold on
end
if ishape==1
    for k=1:ms
        filename = sprintf('object%i',k);
        object = load(filename);
        x = object(:,1);
        y = object(:,2);
        plot(x,y,'k-')
        hold on
    end
end
hold off
axis equal
% axis([0 64 0 64])
xlabel('x')
ylabel('y')

if iplotall==1
    
    figure(2)
    surf(xc(ixcs:kip:ixce),yc(jycs:kip:jyce),uc(jycs:kip:jyce,ixcs:kip:ixce))
    
    figure(3)
    surf(xc(ixcs:kip:ixce),yc(jycs:kip:jyce),vc(jycs:kip:jyce,ixcs:kip:ixce))
    
    figure(4)
    surf(xc(ixcs:kip:ixce),yc(jycs:kip:jyce),p(jycs:kip:jyce,ixcs:kip:ixce))
    
    figure(5)
    surf(xc(ixcs:kip:ixce),yc(jycs:kip:jyce),d(jycs:kip:jyce,ixcs:kip:ixce))
    
end