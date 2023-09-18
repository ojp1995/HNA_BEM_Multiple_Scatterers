% testing to see that the integration for the weights is correct in the
% function v1_weights_mid.m

clear all
addpath('/Users/ojp18/OneDrive - University of Reading/PhD/MATLAB/Polynomial_space_approx/more_robust_attempt')


a = 0;  % start of the interval

b = 2*pi;  % end point of the invertval

k = 10;

x1 = -0.5;

x2 = 0.0001;

x3 = pi;

x4 = 2*pi + 4;

% matlabs approxmations

I_mat_1 = integral(@(t) log(k*abs(x1 - t)), a, b );

I_mat_2 = integral(@(t) log(k*abs(x2 - t)), a, b);

I_mat_3 = integral(@(t) log(k*abs(x3 - t)), a, b );

I_mat_4 = integral(@(t) log(k*abs(x4 - t)), a, b);

I_mat_2a = integral(@(t) log(k*(x2 - t)), a, x2) + integral(@(t) log(k*(t - x2)), x2, b );

I_mat_3a = integral(@(t) log(k*(x3 - t)), a, x3) + integral(@(t) log(k*(t - x3)), x3, b );

% analytically computed answers
I_op_1 = v1_weights_mid(a, b, k, x1);

I_op_2 = v1_weights_mid(a, b, k, x2);

I_op_3 = v1_weights_mid(a, b, k, x3);

I_op_4 = v1_weights_mid(a, b, k, x4);

% comparison

err1 = abs(I_mat_1 - I_op_1);
err2 = abs(I_mat_2 - I_op_2);
err3 = abs(I_mat_3 - I_op_3);
err4 = abs(I_mat_4 - I_op_4);

err2a = abs(I_mat_2a - I_op_2);

err3a = abs(I_mat_3a - I_op_3);

[err1 err2 err3 err4]

[err2a err3a]

disp('Have checked this with wolfram alpha and it would seem that the exact answers given for case 2 and 3 are infact correct')



