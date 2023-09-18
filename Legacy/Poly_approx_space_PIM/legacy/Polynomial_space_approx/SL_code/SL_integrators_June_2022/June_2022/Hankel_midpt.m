function z = Hankel_midpt(a,b,k,s,N)
%Hankel_midpt:
%   Evaluate (1i/4)*\int_a^b H_0^{(1)}(k|s-t|) dt 
%   using simple midpt rule with N intervals,
%   no attempt to deal with singularity in a sophisticated way
h=(b-a)/N;
z=(1i*h/4)*sum(besselh(0,k*abs(s-(a+((1:N)-0.5))*h)));
end

