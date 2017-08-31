%clear
close all

m=8193;
n=8193;
[p]=MPI_IO('p.dat',m,n);
[uc]=MPI_IO('uc.dat',m,n);
[vc]=MPI_IO('vc.dat',m,n);
[rhsp]=MPI_IO('rhsp.dat',m,n);
[ph]=MPI_IO('ph.dat',m,n);
[wo]=MPI_IO('wo.dat',m,n);

[u]=MPI_IO('u.dat',m,n);
[v]=MPI_IO('v.dat',m,n);

[d]=MPI_IO('d.dat',m,n);

wo=wo';
ph=ph';
% figure(1), 
% surf(p) 
% title('p')
% 
% figure(2), 
% surf(uc)
% title('uc') 
% 
% figure(3),
% surf(vc)
% title('vc')
% 
% figure(4),
% surf(rhsp)
% title('rhsp')