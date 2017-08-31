function [matrix,fid,s]=MPI_IO_INT(filename,m,n)

format long
fid=fopen(filename);
% fid=fread(fid,'int16');
matrix=fread(fid,[m n],'int16');
s=size(fid);

% matrix=zeros(m,n);
% for j=1:n
%     for i=1:m
%         index=i+(j-1)*n;
%         matrix(i,j)=fid(index);
%     end
% end
% 
% for i=1:m
%     for j=1:n
%         index=j+(i-1)*m;
%         matrix(i,j)=fid(index);
%     end
% end
% matrix=matrix';

end