% In this script we will be testing a regular midpoint compared to the PIM,
% this time with a general function added into the integral.

clear all

a = 0;  % start of integration
b = 2*pi;  % end of integration

s = 0.4;  % the collocation point we are evaluating the integral at

k = 10;  % wavenumber

C1 = 1; C2 = 2*pi; % constants for smoothing function

N_final  = 20;  % dof = 2^{1:N_final}

f = @(t) exp(1i*k*t);%sin(k*t);

N_true = 25; 
N = 2^N_true;
h = (b - a)/N;
t = [a + h/2:h: b - h/2];  % integration nodes
a_disc = a + (0:N-1)*h;
b_disc = a + (1:N)*h;
I_true_PIM = sum(m1_vec_test(k, s, t, C1, C2).*v1_weights_mid_vec_test(a_disc,b_disc,k,s).*f(t) ...
        + h*m2_vec_test(k, s, t, C1, C2).*f(t));
    
I_mat = 1i*integral(@(y) besselh(0, k*abs(s - y)).*exp(1i*k*y), a, b )/4
for j = 1:N_final
   
    N = 2^j;  % number of dof
    h = (b-a)/N;  %step size
    t = [a + h/2:h: b - h/2];  % integration nodes
    a_disc = a + (0:N-1)*h;
    b_disc = a + (1:N)*h;
    
    I_PIM(j) = sum(m1_vec_test(k, s, t, C1, C2).*v1_weights_mid_vec_test(a_disc,b_disc,k,s).*f(t) ...
        + h*m2_vec_test(k, s, t, C1, C2).*f(t));
    
    I_midpoint(j) = sum(1i*h*f(t).*besselh(0, k*abs(s - t))/4);
    
end

comp_methods = [I_PIM.' I_midpoint.'];

err_PIM = abs(I_true_PIM - I_PIM);

err_mid = abs(I_true_PIM - I_midpoint);

EOC_PIM = log2(err_PIM(1:N_final-1)./err_PIM(2:N_final));
EOC_mid = log2(err_mid(2:N_final)./err_mid(1:N_final-1));

err_comp = [err_PIM.' err_mid.'];

EOC_comp = [EOC_PIM.' EOC_PIM.'];





