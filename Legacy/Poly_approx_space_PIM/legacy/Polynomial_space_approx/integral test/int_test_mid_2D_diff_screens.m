% this function will be adapting the current code so that it will work for
% the screen case, herer we will be integrating midpoint only for x\in
% \Gamma_{1} and y \in \Gamma_{2}.

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

k = 10;  % wavenumber

f0 = @(t) 1;
f1 = @(t) sin(k*t);
f2 = @(t) cos(k*t);
f3 = @(t) exp(1i*k*t);
f4 = @(t) t.^(-1/2) + (L1 - t).^(-1/2);


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

f = f4;

% here we compute the point we are evaluating the integrals at on the
% opposite screen
s = 0.5;  % position along the parameterised screen
s1x = G1(1) + s*(G1(3) - G1(1))/L1;
s1y = G1(2) + s*(G1(4) - G1(2))/L1;

s2x = G2(1) + s*(G2(3) - G2(1))/L2;
s2y = G2(2) + s*(G2(4) - G2(2))/L2;

% "true solution"
I_true_1 = sum(1i*h1*besselh(0, k*sqrt( (s2x - x1 ).^2 + (s2y - y1).^2 )).*f(t1_mid)/4);  % integrating over \Gamma_{1}
I_true_2 = sum(1i*h2*besselh(0, k*sqrt( (s1x - x2 ).^2 + (s1y - y2).^2 )).*f(t2_mid)/4);  % integrating over \Gamma_{2}
%matlab solution
I_mat_1 = 1i*integral(@(t) besselh(0, k*sqrt( (s2x - ( G1(1) + t*(G1(3) - G1(1))/L1 ) ).^2 ...
    + (s2y - ( G1(2) + t*(G1(4) - G1(2))/L1 )).^2 )).*f(t), 0, L1)/4;

I_mat_2 = 1i*integral( @(t) besselh(0, k* sqrt( (s1x - ( G2(1) + t*(G2(3) - G2(1))/L2 )).^2 ...
    + (s1y - (G2(2) + t*(G2(4) - G2(2))/L2) ).^2)).*f(t), 0, L2 )/4;


int_comp_mid_comp_matlab = [I_true_1 I_mat_1 I_true_2 I_mat_2]
%%
for j = 1:N_final
    N = 2^j;
    N1 = N;
    N2 = N;
    [x1, y1, t1, t1_mid, h1, hvector1, L1] = Ninput_discretisation_variables(G1, N1, k);
    [x2, y2, t2, t2_mid, h2, hvector2, L2] = Ninput_discretisation_variables(G2, N2, k);
    
    I_1(j) = sum(1i*h1*besselh(0, k*sqrt( (s2x - x1 ).^2 + (s2y - y1).^2 )).*f(t1_mid)/4);  % integrating over \Gamma_{1}
    I_2(j) = sum(1i*h2*besselh(0, k*sqrt( (s1x - x2 ).^2 + (s1y - y2).^2 )).*f(t2_mid)/4);  % integrating over \Gamma_{2}
 
    I_test_1(j) = midpoint_hankel_f_2D(k, s2x, s2y, x1, y1, t1_mid, h1, f);
    I_test_2(j) = midpoint_hankel_f_2D(k, s1x, s1y, x2, y2, t2_mid, h2, f);
    
end

err_1 = abs(I_true_1 - I_1);
err_2 = abs(I_true_2 - I_2);

EOC_1 = log2(err_1(2:N_final)./err_1(1:N_final-1));
EOC_2 = log2(err_2(2:N_final)./err_2(1:N_final-1));

err_comp = [err_1.' err_2.'];

EOC_comp = [EOC_1.' EOC_2.'];
