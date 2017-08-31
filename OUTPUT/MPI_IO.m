function [matrix]=MPI_IO(filename,m,n)

format long
fid=fopen(filename);
matrix=fread(fid,[m n],'double');

% format long
% fid=fopen(filename);
% fid=fread(fid,'double');
% s=length(fid);
% 
% matrix=zeros(m,n);
% for j=1:n
%     for i=1:m
%         index=i+(j-1)*n;
% % for i=1:m
% %     for j=1:n
% %         index=j+(i-1)*m;
%         matrix(i,j)=fid(index);
%     end
% end

end