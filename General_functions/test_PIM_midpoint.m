% testing the PIM method

clear all

a = 0;
b = 2*pi;
k = 10;

s = pi + 0.1;

Q = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768];

% test f functions
f = @(t) 1; 
f = @(t) sin(k*t);
f = @(t) exp(1i*k*t);

fun = @(t) 1i*besselh(0, k*abs(s - t)).*f(t)/4;

int_mat = integral( @(t) fun(t), a, b);

for q = 1:length(Q)
    h = (b - a)/(Q(q) + 1);
    nq = [a + h/2: h:b - h/2];
    tq = [a:h:b];

    fnq = f(nq);

    I(q) = PIM_int_hankel_f(k, s, h, nq, fnq, tq);

    I_mid(q) = sum(h*fun(nq));

end

% 'True solution'
Q_large = 10*Q(end);
h = (b - a)/(Q_large + 1);
nq = [a + h/2: h:b - h/2];
tq = [a:h:b];

fnq = f(nq);

I_comp = PIM_int_hankel_f(k, s, h, nq, fnq, tq);

err = I_comp - I;

for q = 1:(length(Q) - 1)

    EOC(q) = log2( err(q)/err(q+1));


end

err.'
    
EOC.'

[I.' I_mid.']

int_mat




