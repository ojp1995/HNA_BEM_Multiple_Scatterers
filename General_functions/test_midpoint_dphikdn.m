% testing midpoint_dphikdn_f_diff_screen.m
clear all

a = 0;
b = 2*pi;

k = 10;

s = pi + 0.1;

% Q = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768];
Q_max = 20;
% test f functions
f = @(t) 1; 
f = @(t) sin(k*t);
f = @(t) exp(1i*k*t);

x1 = -1; %[-1, -4, -0.1];
x2 = 1; %[-1, -4, -0.1];

n = [1/sqrt(2), 1/sqrt(2)];

for q = 1:Q_max
% for q = 1:length(Q)
%     h = b/(Q(q) + 1);

    Q = 2^q;
    h = b/(Q + 1);
    nq = [a + h/2: h:b - h/2];

    y1t = nq;
    y2t = nq/2;

    fnq = f(nq);

    I(q) = midpoint_dphikdn_f_diff_screen(k, x1, x2, h, y1t, y2t, fnq, n);

end

Q_large = (2^20)*10;
h = b/(Q_large + 1);
nq = [a + h/2: h:b - h/2];
y1t = nq;
y2t = nq/2;

fnq = f(nq);

I_large_Q = midpoint_dphikdn_f_diff_screen(k, x1, x2, h, y1t, y2t, fnq, n);

err = I_large_Q - I;

% convergence

for q = 1:(length(I) - 1)

    EOC(q) = log2( err(q)/err(q+1));


end

err.'
    
EOC.' 

I.'
