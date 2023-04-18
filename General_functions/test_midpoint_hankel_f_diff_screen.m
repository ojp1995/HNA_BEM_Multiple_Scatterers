% in this script we will be testing that midpoint_hankel_f_diff_screen
% works as expected
clear all
% general variables needed
a = 0;
b = 2*pi;

k = 10;

Q = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024];

x1 = -1%, -4, 0.1];
x2 = 1%, -2, -0.01];


% test f functions
% f = @(t) 1; 
% f = @(t) sin(k*t);
f = @(t) exp(1i*k*t);

% fun = 
% 
% I_mat = 1i*integral(besselh(0, k*( sqrt( (x1 - y1)^2 (x2 - y2)^2 ) ))  )/4



for q = 1:length(Q)
    h = b/(Q(q) + 1);
    nq = [a + h/2: h:b - h/2];

    y1t = nq;
    y2t = nq/2;

    fnq = f(nq);

    I(q) = midpoint_hankel_f_diff_screen(k, x1, x2, y1t, y2t, h, nq, fnq);

    

end
Q_large = Q(end)*10;
h = b/(Q_large + 1);
nq = [a + h/2: h:b - h/2];
y1t = nq;
y2t = nq/2;

fnq = f(nq);


I_large_Q = midpoint_hankel_f_diff_screen(k, x1, x2, y1t, y2t, h, nq, fnq);

err = I_large_Q - I;

% convergence

for q = 1:(length(Q) - 1)

    EOC(q) = log2( err(q)/err(q+1));


end

err.'
    
EOC.' 

I.'


