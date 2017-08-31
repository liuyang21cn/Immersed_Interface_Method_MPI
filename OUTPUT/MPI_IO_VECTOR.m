function [fid,vector,s]=MPI_IO_VECTOR(filename,n)

format long
fid=fopen(filename);
fid=fread(fid,'double');
s=size(fid);

vector=zeros(n,1);
for i=1:n
   vector(i)=fid(i);
end

end