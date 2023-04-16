function I = PIM_int_hankel_f(k, s, h, nq, fnq, tq, C1, C2)
% In this function we will be applying the product integration method for
% the integral:
% \frac{\ri}{4} \int_{0}^{L} H_{0}^{(1)}( k \vert s - t \vert) f(t) \dd t 
%
%
% Inputs:
%
% Outputs:

t_lower = tq(1:end - 1);
t_upper = tq(2:end);

I = 0;

for j = 1:length(s)

    dist = abs(s(j) - nq);

    I = I + sum( -smoothing_function(nq, C1, C2).*besselj(0, k*dist).*fnq.*w1_weights(k, s(j), t_lower, t_upper)/(2*pi) ...
        +  h*fnq.*( 1i*besselh(0, k*dist)/4 + smoothing_function(nq, C1, C2).*besselj(0, k*dist).*log(k*dist)/(2*pi) ) );

   

end