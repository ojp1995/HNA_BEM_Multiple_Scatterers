%Test code, for testing Hankel_midpt_PIM.m and Hankel_midpt.m

clear
a=0;b=2*pi;k=10;s=7*pi/19;C1=1;C2=pi; %test parameters
j=25; N=2^j; %to compute exact solution
x=Hankel_midpt(a,b,k,s,N);y=Hankel_midpt_PIM(a,b,k,s,C1,C2,N);
xexact=x;yexact=y; % "exact" solution = 1.2369e-03 + 4.4904e-02i for midpt and PIM
M=20; xapprox=zeros(M+1,1);yapprox=zeros(M+1,1); %setting up loop for approx solns
for j=0:M;
  N=2^j;
  x=Hankel_midpt(a,b,k,s,N);y=Hankel_midpt_PIM(a,b,k,s,C1,C2,N);
  xapprox(j+1)=x;yapprox(j+1)=y;
end
xerror=abs(yexact-xapprox);yerror=abs(yexact-yapprox); %compute errors - PIM as "exact"
disp([xerror yerror]) %display errors:  midpt ~ 10^{-7}, PIM ~ 10^{-12}
log2(xerror(2:21)./xerror(1:20)) %EOC for midpt variable
log2(yerror(2:21)./yerror(1:20)) %EOC for PIM = -2.0 => quadratic convergence

