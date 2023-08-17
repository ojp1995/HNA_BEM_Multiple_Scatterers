function I = graded_rescalled_PIM_int_hankel_f(k, s, w, nq, fnq, tq, C1,...
    C2, L)
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

I = zeros(length(s), 1);
for j = 1:length(s)

%     dist = abs(s(j) - nq);

    I(j, 1) = sum((m1(k, s(j), L - nq, C1, C2).*w1_weights_rescalled(k, ...
        s(j), t_lower, t_upper, L) + w.*m2(k, s(j), L - nq, C1, C2)).*fnq );
end