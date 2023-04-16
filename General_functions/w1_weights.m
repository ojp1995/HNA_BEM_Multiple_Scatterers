function w = w1_weights(k, s, tj_1, tj)
% In this function we will be computing the integral which are the
% analytical weights for the PIM method for integrating the hankel function
% as a kernel.
%
%
% Inputs:
%
% Outputs:

w = (s - tj_1).*log( k*abs(s - tj_1) ) + (tj - s).*log(k*abs( tj - s )) + tj_1 - tj;
