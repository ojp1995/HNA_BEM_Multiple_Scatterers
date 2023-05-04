function z = m2(k, s, nq, C1, C2)
% In this function we will be computing the following kernel:
% m2(s, nq) = i H_{0}^{(1)}(k | s - nq |) - \sigma_{1}(s, t)m_{1}(s, t)

% select = (s ~= nq);

select = (s == nq);

dist = s - nq;

z = 1i*besselh(0, k*dist)/4 ...
    + smoothing_function(nq, C1, C2).*besselj(0, k*dist).*log(k*dist)/(2*pi); 

z(select) = 1i/4 - 2*log(1/2)/pi;

