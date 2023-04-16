% testing the product integration method for a range of functions.

clear all

a = 0;  % start of interval
b = 1;  % end of interval

N = [5, 10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240];  % number of intervals

k = 100;

m = @(x) sin(x); 

% Test 1
phi = @(x) 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TEST 1a x > a %%%%%%%%%%%%%%%%
x1 = -0.5;
I_true1 = -0.011165050397663045277015931956078230060779551951865292685125177;

% I_mid1 = phi(x1)*v1_weights_mid(a, b, x1);
for j = 1:length(N)
    I_PIM1(j) = ProductMidpoint_log(a, b, x1, m, phi, k, N(j));
    
    err1(j) =  abs(I_PIM1(j) - I_true1);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% TEST 1b x \in(a, b) %%%%%%%%%%
x2 = 0.3;
I_true2 = 0.127711;
% I_mid2 = phi(x2)*v1_weights_mid(a, b, x2);
for j= 1:length(N)
    I_PIM2(j) = ProductMidpoint_log(a, b, x2, m, phi, k, N(j));
    
    err2(j) = abs(I_PIM2(j) - I_true2);
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TEST 1c x > b %%%%%%%%%%%%%%%%
x3 = 1.4;
I_true3 = -0.060332945499480093966004087240104760059390470340480776583700926;
%I_mid3 = phi(x3)*v1_weights_mid(a, b, x3);
for j = 1:length(N)
    I_PIM3(j) = ProductMidpoint_log(a, b, x3, m, phi, k, N(j));
    
    err3(j) = abs(I_PIM3(j) - I_true3);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TEST 1d x >> b %%%%%%%%%%%%%%%%
x4 = 1e6;
I_true4 = -10.0181;
% I_mid4 = phi(x4)*v1_weights_mid(a, b, x4);
for j = 1:length(N)
    I_PIM4(j) = ProductMidpoint_log(a, b, x4, m, phi, k, N(j));
    
    err4(j) =  abs(I_PIM4(j)  - I_true4);
end



% computing EOC (I think)

for j = 1:(length(N) - 1)
    
    EOC1(j) = log2( err1(j)/err1(j+1) );
    
    EOC2(j) = log2( err2(j)/err2(j+1) );
    
    EOC3(j) = log2( err3(j)/err3(j+1) );
    
    EOC4(j) = log2( err4(j)/err4(j+1) );

end

T_EOC = [EOC1.' EOC2.' EOC3.' EOC4.'];

T_err = [err1.' err2.' err3.' err4.'];