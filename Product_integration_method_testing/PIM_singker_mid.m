%PIM using midpoint rule test, singular kernel
clear all
a = 0;  % start of interval
b = 1;

lnz_true = -1/4;

N = [5, 10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240];  % number of intervals

% function we are integrating
sigma1 = @(x) log(x);
m1 = @(x) x; %x.^2;

sigma2 = @(x) 1;
m2 = @(x) 0;

for j = 1: length(N)
    
    h = (b-a)/N(j);  % step size
    
    x = [a: h: b];  % discretisation of [a, b]
    
    %weights computation
    w(1) = h*log(h) - h;
    for l = 2:N(j)
        w(l) = x(l+1)*(log(x(l+1)) - 1 ) - x(l)*(log(x(l)) - 1 );
    end

    n = [a+h/2:h: b - h/2];  % these are the midpoints of each interval, the nodes
    
    z(j) = PIM_midpoint(sigma1, m1, sigma2, m2, w, n);
%     disp(w.')
%     clear w x
    
    h_square_inv(j) = 1/(h^2);

    
end

disp(z.')
% disp(lnz_true)

err = z - lnz_true;

[err.' h_square_inv.']



% need to see if the error behaves like h^{-2}
for j = 1:(length(N) - 1)
    
    EOC(j) = log2(err(j)/err(j+1));
    
end


disp(EOC.')
