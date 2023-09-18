% in this function we will look at whether the integration method we are
% using actually works using the product integration method.
%
% We are going to evlauate the integral 
% \int_{a}^{b} H_{0}^{(1)}(k |x - t|)dt, for some t \in [a, b]. This will
% be tested for finer and finer mesh until the integral stabilises.

clear all

addpath('/Users/ojp18/OneDrive - University of Reading/PhD/MATLAB/Polynomial_space_approx/more_robust_attempt')

a = 0; % begining of the integal
b = 2*pi; % end point of the integral

k = 1;  % hardwired in at some points

L = b - a;

x = 2;  % point we have the singularity at

C1 = 1;
C2 = pi;

% matlab solve
I_mat = 1i*integral( @(t) besselh(0, abs(x - t)), a, b )/4;

I_mat_sing = -integral(@(t) smoothing_fun(t, C1, C2).*log(abs(x - t)).*besselj(0, abs(x - t)), a, b)/(2*pi);

I_mat_smooth = integral( @(t) 1i*besselh(0, abs(x - t))/4 +  smoothing_fun(t, C1, C2).*log(abs(x - t)).*besselj(0, abs(x - t))/(2*pi), a, b);

I_mat_total = I_mat_sing + I_mat_smooth;

N_final = 15;  % discretisation of screen will be 2^[1:1:N_final]



I = zeros(N_final, 1);
I_sing_total = zeros(N_final, 1);
I_smooth_total = zeros(N_final, 1);
I_besselj = zeros(N_final, 1);
for j = 1:N_final % looping through the different discretisations
    
    N = 2^j;  % discretisation parameter
    
    h = L/N;  % step size
    
    t = [a:h:b];  % discretisation of screen for integration
    
    n = [a+h/2: h: b - h/2]; % midpoints of the interval, our integration nodes
    
    % evaluating the integral
    for l = 1:(length(t)-1)
        
        dist = abs( x - n(l) );
%         I_smoothfun_test(j, 1) = I_smoothfun_test(j, 1) +
%         smoothing_fun(dist, C1, C2); % test not needed, it will grow as
%         the number of intervals are increased.

        

%         I_besselj(j, 1) = I_besselj(j, 1) + besselj(0, dist);

           
        
        I_sing = m1_tilde(1, dist, C1, C2)*v1_weights_mid(t(l), t(l+1), 1, n(l));
        
        I_sing_total(j, 1) = I_sing_total(j, 1) + I_sing;
        
        I_smooth = h*m2(1, dist, C1, C2);
        I_smooth_total(j, 1) = I_smooth_total(j, 1) + I_smooth;
        
        I(j, 1) = I(j, 1) + I_sing + I_smooth;
        
    end
    
    I_new_smooth(j) = 1i*h*sum( besselh(0, k*abs(x - ( a + ([1:N]-1/2)*h ) )))/4 ...
        + sum(smoothing_fun( (a + ([1:N] - 1/2)*h) , C1, C2).*log(k*abs( x - ( a + ([1:N]-1/2)*h ) )).*besselj(0,k*abs( x - ( a + ([1:N]-1/2)*h )))  )/(2*pi); % still need to add in smoothing function and log and J_0 
    % comparison to matlab solver
    
    err(j) = abs(I_mat - I(j, 1))/abs(I_mat);
    
    
    
    
end

[I_mat I.'].'
