% b=load('jcp_rhs.txt');
% bb = b(1:256);
% [ sol,alfa,beta,r1,w1 ] = ConjugateGradient(bb);

load jcp1.dat
bb = jcp1(1:256);
[ sol,alfa,beta,r1,w1,rnorm,w0,r0] = ConjugateGradient(bb);
