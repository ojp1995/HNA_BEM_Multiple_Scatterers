function X = smoothing_function(t, C1, C2)
% This function is multiplied by the smooth function to smoothly reduce it
% to zero so that it equals 1 when |t|=< 1, 0 when C2 <= |t| and equal to
% the 1 + exp( (pi - 1)/(pi - |t|) - (1 - pi)/(1 - |t|) )^{-1} when 
% C1 < |t| < C2
% 
%
% Problem parmeters:
% t is a vector of points we are inputting to the function
% C1 is the first constant, advised to be 1
% C2 is the second constant, advised to be pi


X = zeros(length(t), 1);

% truth table method
a = abs(t)<=C1; % vector of 1s for the case |t|<= C1
b = abs(t)>=C2; % vector of 1s for the case C2<= |t|
c = 1 - a - b; % inbetween case

% here computing vector, c*inbetween case + a*C1 case, b not needed as only
% 0's remain.
X = c.*1./( 1 + exp( (C2 - C1)./( C2 - abs(t)) - (C1 - C2)./(C1 - abs(t))) ) + a;
