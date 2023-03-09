% In this script we will be testing that the Product integration methods
% works as expected for singular and oscailltory integrals

% first testing the weights are performing correctly:

clear all
a = 0;  % lower bound of integral
b = 1;  % upper bound of integral
k = 1;  % k not needed yet but may be in future

%% Test 1 x > a
x1 = - 0.5;
I_true1 = -0.0452287;

I_approx1 = v1_weights_mid(a,b,k,x1);

err1 = I_true1 - I_approx1;

%% Test 2 x = a
x2 = a;
I_true2 = -1;
I_approx2 = v1_weights_mid(a,b,k,x2);

err2 = I_true2 - I_approx2;


%% TEST 3 x \in(a, b) %%%%%%%%%%
x3 = 0.3;
I_true3 = - 1 + log( ( 3^(3/10) )*( 7^(7/10) )/10);
I_approx3 = v1_weights_mid(a, b, k, x3);

err3 = I_true3 - I_approx3;


%% TEST 4 x = b %%%%%%%%%%%%%%%%
x4 = b;
I_true4 = -1;
I_approx4 = v1_weights_mid(a, b, k, x4);

err4 = I_true4- I_approx4;



%% TEST 5 x > b %%%%%%%%%%%%%%%%
x5 = 1.4;
I_true5 = -0.162423;
I_approx5 = v1_weights_mid(a, b, k, x5);

err5 =  I_true5- I_approx5;


%% TEST 6 x >> b %%%%%%%%%%%%%%%%
x6 = 1e6;
I_true6 = 13.816;

I_approx6 = v1_weights_mid(a, b, k, x6);

err6 =  I_true6- I_approx6;

err = [err1, err2, err3, err4, err5, err6];

disp(err)



