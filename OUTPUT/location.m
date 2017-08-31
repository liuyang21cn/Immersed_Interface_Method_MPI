function location(m,n,rhs1,rhs2)

max=-10000000000000000000000;
for i=1:m
    for j=1:n
        err=abs(rhs1(i,j)-rhs2(i,j));
        if err > max
            imax = i;
            jmax = j;
            max = err;
        end
    end
end

fprintf('max = %.8f, imax= %i, jmax= %i \n', max, imax, jmax);

end