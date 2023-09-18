function z = Hankel_midpt_PIM(a,b,k,s,C1,C2,N)
%Hankel_midpt_PIM:
%   Evaluate (1i/4)*\int_a^b H_0^{(1)}(k|s-t|) dt 
%   using product integration midpt rule with N intervals,

h=(b-a)/N;
nodes=a+((1:N)-0.5)*h;
z=sum((v1_weights_mid(a+(0:N-1)*h,a+(1:N)*h,k,s).*...
    m1_tilde(k,s,nodes,C1,C2))+h*m2(k,s,nodes,C1,C2));
end

