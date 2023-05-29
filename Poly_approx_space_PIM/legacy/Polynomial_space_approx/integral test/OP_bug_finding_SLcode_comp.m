
% does code work??? _ BUG HUNTING!!
%
% In this piece of code we will be looking to try and find the bugs in the
% code that I have written and will be comparing it to Steve's good code

clear all

a=0;b=2*pi;k=10;x=7*pi/19;C1=1;C2=pi; %test parameters
N = 2^15;

addpath('/Users/ojp18/OneDrive - University of Reading/PhD/MATLAB/Polynomial_space_approx/more_robust_attempt')

L = b-a;
h = L/N;  % step size

t = [a:h:b];  % discretisation of screen for integration

n = [a+h/2: h: b - h/2]; % midpoints of the interval, our integration nodes

N_final = 1;
I = zeros(N_final, 1);
I_sing_total = zeros(N_final, 1);
I_smooth_total = zeros(N_final, 1);
I_besselj = zeros(N_final, 1);
m1_tilde_ans = 0;
m2_tilde_ans = 0;
weights = 0;
% evaluating the integral
for j = 1:1

    for l = 1:(length(t)-1)

        dist = abs( x - n(l) );

        m1_tilde_ans = m1_tilde_ans + m1_tilde(k, dist, C1, C2);
        
        weights = weights +v1_weights_mid(t(l), t(l+1), k, n(l));
        
        m2_tilde_ans = m2_tilde_ans + m2(k, dist, C1, C2);

        I_sing = m1_tilde(k, dist, C1, C2)*v1_weights_mid(t(l), t(l+1), k, n(l));

        I_sing_total(j, 1) = I_sing_total(j, 1) + I_sing;

        I_smooth = h*m2(k, dist, C1, C2);
        I_smooth_total(j, 1) = I_smooth_total(j, 1) + I_smooth;

    end
end

m1_tilde_ans
m2_tilde_ans
weights