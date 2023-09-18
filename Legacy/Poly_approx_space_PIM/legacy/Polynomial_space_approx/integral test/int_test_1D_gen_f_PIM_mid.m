% In this script we will be testing a regular midpoint compared to the PIM,
% this time with a general function added into the integral.

clear all

a = 0;  % start of integration
b = 2*pi;  % end of integration

s = pi+0.15;  % the collocation point we are evaluating the integral at

k = 10;  % wavenumber

C1 = 1; C2 = 2*pi; % constants for smoothing function

N_final  = 20;  % dof = 2^{1:N_final}

f = @(t)  exp(1i*k*t); % cos(k*t);

N_true = 25; 
N = 2^N_true;
h = (b - a)/N;
t_mid = [a + h/2:h: b - h/2];  % integration nodes
a_disc = a + (0:N-1)*h;
b_disc = a + (1:N)*h;
I_true_PIM = sum(m1_vec_test_genf(k, s, t_mid, C1, C2, f).*v1_weights_mid_vec_test(a_disc,b_disc,k,s) ...
        + h*m2_vec_test_genf(k, s, t_mid, C1, C2, f));

I_mat = 1i*integral(@(y) besselh(0, k*abs(s - y)).*f(y), a, b )/4
    
I_mat1 = 1i*integral(@(y) besselh(0, k*abs(s - y)), a, b )/4;
I_mat2 = 1i*integral(@(y) besselh(0, k*abs(s - y)).*sin(k*y), a, b )/4;
I_mat3 = 1i*integral(@(y) besselh(0, k*abs(s - y)).*cos(k*y), a, b )/4;
I_mat4 = 1i*integral(@(y) besselh(0, k*abs(s - y)).*exp(1i*k*y), a, b )/4;

%%
for j = 1:N_final
   
    N = 2^j;  % number of dof
    h = (b-a)/N;  %step size
    t = [a:h:b];
    t_mid = [a + h/2:h: b - h/2];  % integration nodes
    a_disc = a + (0:N-1)*h;
    b_disc = a + (1:N)*h;
    
    I_PIM_test(j) = PIM_mid_hankel(k, s, t_mid, t, C1, C2, f, h);
    
    I_PIM(j) = sum(m1_vec_test_genf(k, s, t_mid, C1, C2, f).*v1_weights_mid_vec_test(a_disc,b_disc,k,s) ...
        + h*m2_vec_test_genf(k, s, t_mid, C1, C2, f));
    
    I_midpoint(j) = sum(1i*h*f(t_mid).*besselh(0, k*abs(s - t_mid))/4);
    
end

comp_methods_PIM_mid = [I_PIM_test.' I_midpoint.']
% comp_methods = [I_PIM_test.' I_PIM I_midpoint.'] % old, I_PIM and
% I_PIM_test are the same

err_PIM = abs(I_true_PIM - I_PIM);

err_mid = abs(I_true_PIM - I_midpoint);

EOC_PIM = log2(err_PIM(2:N_final)./err_PIM(1:N_final-1));
EOC_mid = log2(err_mid(2:N_final)./err_mid(1:N_final-1));

err_comp_PIM_mid = [err_PIM.' err_mid.'];

EOC_comp_PIM_mid = [EOC_PIM.' EOC_mid.'];





