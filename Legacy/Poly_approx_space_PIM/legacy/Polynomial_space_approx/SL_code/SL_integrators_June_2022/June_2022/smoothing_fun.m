function X = smoothing_fun(t, C1, C2)
% Smoothing function:
% t is a vector of points we are inputting to the function
% C1 is the first constant, advised to be 1
% C2 is the second constant, advised to be pi
a = abs(t)<=C1; % vector of 1s for the case |t|<= C1
b = abs(t)>=C2; % vector of 1s for the case C2<= |t|
c = 1 - a - b; % c=1 if C1<t<C2
X = c.*1./(1+exp((C2-C1)./(C2-abs(t))-(C1-C2)./(C1-abs(t)))) + a;
end