function [ sol,alfa,beta,r1,w1,rnorm,w_old,r_old ] = ConjugateGradient(b)

n=length(b);
A=diag(-2*ones(n,1),0)+diag(ones(n-1,1),1)+diag(ones(n-1,1),-1);
A(1,n)=1;
A(n,1)=1;

sol=zeros(n,1);
r0=b-A*sol;
r1=zeros(n,1);
w1=zeros(n,1);
w0=r0;
rnorm = zeros(100,1);

for k=1:31
    
    w_old = w0;
    r_old = r0;
    alfa=(r0'*r0)/(w0'*A*w0);
    if k==32
        alfa
        r0'*r0
        w0'*A*w0
    end
%     r'*r
%     w'*A*w
    sol = sol+alfa*w0;
    r1 = r0-alfa*A*w0;
    rnorm(k) = norm(r1,2);
    norm(r1,2);
    if norm(r1,2)<1e-10
        k
        norm(r1,2)
        break
    else
        beta = (r1'*r1)/(r0'*r0);
        if k==31
            beta
        end
        w1=r1+beta*w0;
    end
    w0=w1;
    r0=r1;

end
