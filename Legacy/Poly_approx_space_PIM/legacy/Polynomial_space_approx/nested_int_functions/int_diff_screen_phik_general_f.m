function I = int_diff_screen_phik_general_f(x_int, y_int, x_col, y_col, k, FUN, h)
% This function will integrate using the midpoint rule of the kernel
% \phik(x, y) f(y) dy. 
% 
%
% Problem parameters
% (x1, y1) is the coordinates of G1
% (x2, y2) are the coordinates of G2
% k is the wavenumber
%
% Discretisation parameters:
% h are the integration weights
%

dist = sqrt( (y_col - y_int)^2 + (x_col - x_int)^2 );

I = 1i*h*besselh(0, k*dist )*FUN/4;