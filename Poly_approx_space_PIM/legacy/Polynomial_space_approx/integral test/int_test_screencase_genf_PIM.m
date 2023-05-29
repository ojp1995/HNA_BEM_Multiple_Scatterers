% this function will be adapting the current code so that it will work for
% the screen case 

clear all


% introduction given for screen currently then we will look at solving for
% a specific configuration

% discretisation 
N_true = 25;  % true number of dof 2^{N_true}
N_final = 20;  % number of dof 2^{1:N_final}

G1 = [-2*pi, 2*pi, 0, 0];
L1 = sqrt( (G1(3) - G1(1))^2 + (G1(4) - G1(2))^2 );  % length of G1

G2 = [ 2*pi, 0, 5*pi, 3*pi];
L2 = sqrt( (G2(3) - G2(1))^2 + (G2(4) - G2(2))^2 );  % length of G2

% C_wl= 1/40

k = 5;  % wavenumber

f0 = @(t) 1;
f1 = @(t) sin(k*t);
f2 = @(t) exp(1i*k*t);


theta = 0;

% constants needed for the smoothing function
C1 = 1;
C2 = pi;

% N1 = ceil(k*L1./(C_wl*2*pi)) % number of itervals on G1
% N2 = ceil(k*L2./(C_wl*2*pi)) % number of intervals on G2
s = 0.5;

N = 2^N_true;
N1 = N;
N2 = N;
[x1, y1, t1, t1_mid, h1, hvector1, L1] = Ninput_discretisation_variables(G1, N1, k);
[x2, y2, t2, t2_mid, h2, hvector2, L2] = Ninput_discretisation_variables(G2, N2, k);
% now we start the solve:

f = f0;

I_PIM_G1_true = PIM_mid_hankel(t1(1), k, s, t1_disc, C1, C2, f, h1, N1);
I_PIM_G2_true = PIM_mid_hankel(t2(1), k, s, t2_disc, C1, C2, f, h2, N2);
%%

% lower order approximations.
for j = 1:N_final
    N = 2^j;
    N1 = N;
    N2 = N;
    
    [x1, y1, t1, t1_mid, h1, hvector1, L1] = Ninput_discretisation_variables(G1, N1, k);
    [x2, y2, t2, t2_mid, h2, hvector2, L2] = Ninput_discretisation_variables(G2, N2, k);
    
    % now we start the solve:
    
    I_PIM_G1(j) = PIM_mid_hankel(t1(1), k, s, t1_disc, C1, C2, f, h1, N1);
    I_PIM_G2(j) = PIM_mid_hankel(t2(1), k, s, t2_disc, C1, C2, f, h2, N2);
    
    I_mid_G1(j) = sum(1i*h1*f(t1_disc).*besselh(0, k*abs(s - t1_disc))/4);
    I_mid_G2(j) = sum(1i*h2*f(t2_disc).*besselh(0, k*abs(s - t2_disc))/4);
    
    
end

comp_methods = [I_PIM_G1.'  I_mid_G1.' I_PIM_G2.' I_mid_G2.']

err_PIM_G1 = abs(I_PIM_G1 - I_PIM_G1_true);
err_PIM_G2 = abs(I_PIM_G2 - I_PIM_G2_true);


I_G2_mat = 1i*integral(@(t)besselh(0, k*abs(s - t)).*f(t), 0, L2)/4 
err_mid_G1 = abs(I_mid_G1 - I_PIM_G1_true);
err_mid_G2 = abs(I_mid_G2 - I_PIM_G2_true);

EOC_PIM_G1 = log2(err_PIM_G1(2:N_final)./err_PIM_G1(1:N_final-1));
EOC_PIM_G2 = log2(err_PIM_G2(2:N_final)./err_PIM_G2(1:N_final-1));

EOC_mid_G1 = log2(err_mid_G1(2:N_final)./err_mid_G1(1:N_final-1));
EOC_mid_G2 = log2(err_mid_G2(2:N_final)./err_mid_G2(1:N_final-1));

err_comp_G1 = [err_PIM_G1.' err_mid_G1.'];
err_comp_G2 = [err_PIM_G2.' err_mid_G2.'];


EOC_comp_G1 = [EOC_PIM_G1.' EOC_mid_G1.'];
EOC_comp_G2 = [EOC_PIM_G2.' EOC_mid_G2.'];
    


