function z = m2(k, s, nq, C1, C2)
% In this function we will be computing the following kernel:
% m2(s, nq) = i H_0^(1)(k |s - nq|)/4 - \sigma_1(s, t)m_{1}(s, t)

% select = (s ~= nq);

select = (s == nq);

dist = s - nq;

z = 1i*besselh(0, k*dist)/4 ...
    + m1(k, s, nq, C1, C2).*log(k*dist)/(2*pi); %% should this be w1_weights??

z(select) = 1i/4 - 2*(log(1/2) - (-psi(1)))/pi;

