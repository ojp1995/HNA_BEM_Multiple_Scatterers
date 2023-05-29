% this is a more explicit test to see if code differes depending on the way 
% the coordinates are input in the method 
%
% SIMONS CODE
%
% We will run the code 4 times, once for the coordinates oriented as in
% other examples, once with the coordinates flipped for G1, another time
% with flipped for G2 and finally with both flipped.
% (By flipped I mean the start coordianted put in the end and the end put 
% in the start etc)

clear all

addpath('/Users/ojp18/OneDrive - University of Reading/PhD/MATLAB/Polynomial_space_approx/Simon_code')

r1start = [-2*pi 2*pi];
r1end = [0 0];

r2start = [4*pi 0];
r2end = [7*pi 3*pi];

L1 = sqrt( (r1end(2) - r1start(2))^2 + ( r1end(1) - r1start(1) )^2 )
L2 = sqrt( (r2end(2) - r2start(2))^2 + ( r2end(1) - r2start(1) )^2 )


r1startA = r1end;
r1endA = r1start;
L1A = sqrt( (r1endA(2) - r1startA(2))^2 + ( r1endA(1) - r1startA(1) )^2 )


r2startB = r2end;
r2endB = r2start;
L2B = sqrt( (r2endB(2) - r2startB(2))^2 + ( r2endB(1) - r2startB(1) )^2 )

C_wl = 1/40;

d = [0,-1]; % the direction of the incident wave

format compact
k = 5 % the wavenumber
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  CONTROL CASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% control case, coordinates oriented as expected:
N1 = ceil(k*L1./(C_wl*2*pi)); N2 = ceil(k*L2./(C_wl*2*pi));% the number of boundary elements on screens 1 and 2

h1 = norm(r1end-r1start)/N1; % the length of each element on screen 1
h2 = norm(r2end-r2start)/N2; % the length of each element on screen 2
% Next calculate the x and y coordinates of the element midpoints on screen
% 1 and then on screen 2

x1 = r1start(1)+((1:N1)-0.5)*(r1end(1)-r1start(1))/N1;
y1 = r1start(2)+((1:N1)-0.5)*(r1end(2)-r1start(2))/N1;
x2 = r2start(1)+((1:N2)-0.5)*(r2end(1)-r2start(1))/N2;
y2 = r2start(2)+((1:N2)-0.5)*(r2end(2)-r2start(2))/N2;
h1vector = h1*ones(size(x1)); % Lengths of the elements on screen 1 
h2vector = h2*ones(size(x2)); % and on screen 2
h = [h1vector,h2vector]; % the lengths of the elements
x = [x1,x2]; % the x-coords of the element midpoints
y = [y1,y2]; % and their y-coords

phi1_con = bem2(x1,y1,h1vector,k,d);
phi2_con = zeros(size(x2));
phi_con = [phi1_con.',phi2_con];

phi10_con = phi1_con; %Keeping the zero iteration

phi2_con = bem2BS(x1,y1,x2,y2,h1vector,h2vector,phi1_con,k,d); % approximations to the values of phi at the element midppoints, phi(j) the value at (x(j),y(j))
phi_con = [phi1_con.',phi2_con.'];

%%% First iteration save
phi21_con = phi2_con; % keeping the 1 iteration


% Now let's compute phi12, the 2nd iteration
phi1_con = bem2BS(x2,y2,x1,y1,h2vector,h1vector,phi2_con,k,d);
phi_con = [phi1_con.',phi2_con.'];

phi12_con = phi1_con; 

% adapting from here on for more iterations

% Now let's compute phi23, the 3rd iteration
phi2_con = bem2BS(x1,y1,x2,y2,h1vector,h2vector,phi1_con,k,d); % approximations to the values of phi at the element midppoints, phi(j) the value at (x(j),y(j))
phi_con = [phi1_con.',phi2_con.'];

phi23_con = phi2_con;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% CASE A: G1 coordinates swapped %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N1 = ceil(k*L1./(C_wl*2*pi)); N2 = ceil(k*L2./(C_wl*2*pi));% the number of boundary elements on screens 1 and 2

h1 = norm(r1endA-r1startA)/N1; % the length of each element on screen 1
h2 = norm(r2end-r2start)/N2; % the length of each element on screen 2
% Next calculate the x and y coordinates of the element midpoints on screen
% 1 and then on screen 2

x1 = r1startA(1)+((1:N1)-0.5)*(r1endA(1)-r1startA(1))/N1;
y1 = r1startA(2)+((1:N1)-0.5)*(r1endA(2)-r1startA(2))/N1;
x2 = r2start(1)+((1:N2)-0.5)*(r2end(1)-r2start(1))/N2;
y2 = r2start(2)+((1:N2)-0.5)*(r2end(2)-r2start(2))/N2;
h1vector = h1*ones(size(x1)); % Lengths of the elements on screen 1 
h2vector = h2*ones(size(x2)); % and on screen 2
h = [h1vector,h2vector]; % the lengths of the elements
x = [x1,x2]; % the x-coords of the element midpoints
y = [y1,y2]; % and their y-coords

phi1_A = bem2(x1,y1,h1vector,k,d);
phi2_A = zeros(size(x2));
phi_A = [phi1_A.',phi2_A];

phi10_A = phi1_A; %Keeping the zero iteration

phi2_A = bem2BS(x1,y1,x2,y2,h1vector,h2vector,phi1_A,k,d); % approximations to the values of phi at the element midppoints, phi(j) the value at (x(j),y(j))
phi_A = [phi1_A.',phi2_A.'];

%%% First iteration save
phi21_A = phi2_A; % keeping the 1 iteration


% Now let's compute phi12, the 2nd iteration
phi1_A = bem2BS(x2,y2,x1,y1,h2vector,h1vector,phi2_A,k,d);
phi_A= [phi1_A.',phi2_A.'];

phi12_A = phi1_A; 

% adapting from here on for more iterations

% Now let's compute phi23, the 3rd iteration
phi2_A = bem2BS(x1,y1,x2,y2,h1vector,h2vector,phi1_A,k,d); % approximations to the values of phi at the element midppoints, phi(j) the value at (x(j),y(j))
phi_A = [phi1_A.',phi2_A.'];

phi23_A = phi2_A;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%Comparison to control case%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
plot(x1, real(phi10_A - phi10_con), 'DisplayName', 'r = 0')
hold on
plot(x1, real(phi10_A), 'DisplayName', 'Case A')
plot(x1, real(phi10_con), 'DisplayName', 'Case control')
plot(x1, flipud(real(phi10_A)), 'DisplayName', 'Case A flipped', 'LineStyle', '-.')
plot(x1, real( phi10_con - flipud(phi10_A) ), 'DisplayName', 'Difference, Case A flipped', 'LineStyle', ':')
title('difference between case A and control case, $\phi_{1}^{(0)}$')
xlabel('x')
ylabel('$\phi_{1}^{(0)}$')
legend show

figure()
plot(x2, real(phi21_A - phi21_con), 'DisplayName', 'r = 1')
title('difference between case A and control case, $\phi_{2}^{(1)}$ ')
xlabel('x')
ylabel('$\phi_{2}^{(1)}$')

figure()
plot(x1, real(phi12_A - phi12_con), 'DisplayName', 'r = 2')
hold on
plot(x1, real(phi12_A), 'DisplayName', 'Case A')
plot(x1, real(phi12_con), 'DisplayName', 'Case control')
plot(x1, flipud(real(phi12_A)), 'DisplayName', 'Case A flipped', 'LineStyle', '-.')
plot(x1, real( phi12_con - flipud(phi12_A) ), 'DisplayName', 'Difference, Case A flipped', 'LineStyle', ':')
title('difference between case A and control case, $\phi_{1}^{(2)}$')
xlabel('x')
ylabel('$\phi_{1}^{(2)}$')
legend show

figure()
plot(x2, real(phi23_A - phi23_con), 'DisplayName', 'r = 3')
title('difference between case A and control case, $\phi_{2}^{(3)}$')
xlabel('x')
ylabel('$\phi_{2}^{(3)}$')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% CASE B: G2 coordinates swapped %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N1 = ceil(k*L1./(C_wl*2*pi)); N2 = ceil(k*L2./(C_wl*2*pi));% the number of boundary elements on screens 1 and 2

h1 = norm(r1end-r1start)/N1; % the length of each element on screen 1
h2 = norm(r2endB-r2startB)/N2; % the length of each element on screen 2
% Next calculate the x and y coordinates of the element midpoints on screen
% 1 and then on screen 2

x1 = r1start(1)+((1:N1)-0.5)*(r1end(1)-r1start(1))/N1;
y1 = r1start(2)+((1:N1)-0.5)*(r1end(2)-r1start(2))/N1;
x2 = r2startB(1)+((1:N2)-0.5)*(r2endB(1)-r2startB(1))/N2;
y2 = r2startB(2)+((1:N2)-0.5)*(r2endB(2)-r2startB(2))/N2;
h1vector = h1*ones(size(x1)); % Lengths of the elements on screen 1 
h2vector = h2*ones(size(x2)); % and on screen 2
h = [h1vector,h2vector]; % the lengths of the elements
x = [x1,x2]; % the x-coords of the element midpoints
y = [y1,y2]; % and their y-coords

phi1_B = bem2(x1,y1,h1vector,k,d);
phi2_B = zeros(size(x2));
phi_B = [phi1_B.',phi2_B];

phi10_B = phi1_B; %Keeping the zero iteration

phi2_B = bem2BS(x1,y1,x2,y2,h1vector,h2vector,phi1_B,k,d); % approximations to the values of phi at the element midppoints, phi(j) the value at (x(j),y(j))
phi_B = [phi1_B.',phi2_B.'];

%%% First iteration save
phi21_B = phi2_B; % keeping the 1 iteration


% Now let's compute phi12, the 2nd iteration
phi1_B = bem2BS(x2,y2,x1,y1,h2vector,h1vector,phi2_B,k,d);
phi_B = [phi1_B.',phi2_B.'];

phi12_B = phi1_B; 

% adapting from here on for more iterations

% Now let's compute phi23, the 3rd iteration
phi2_B = bem2BS(x1,y1,x2,y2,h1vector,h2vector,phi1_B,k,d); % approximations to the values of phi at the element midppoints, phi(j) the value at (x(j),y(j))
phi_B = [phi1_B.',phi2_B.'];

phi23_B = phi2_B;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%Comparison to control case%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
plot(x1, real(phi10_B - phi10_con), 'DisplayName', 'r = 0')
title('difference between case B and control case')
xlabel('x')
ylabel('$\phi_{1}^{(0)}$')

figure()
plot(x2, real(phi21_B - phi21_con), 'DisplayName', 'r = 1')
hold on
title('difference between case B and control case')
plot(x2, real(phi21_B), 'DisplayName', 'Case B')
plot(x2, real(phi21_con), 'DisplayName', 'Case control')
plot(x2, flipud(real(phi21_B)), 'DisplayName', 'Case B flipped', 'LineStyle', '-.')
plot(x2, real( phi21_con - flipud(real(phi21_B))), 'DisplayName', 'Diff with Case B flipped', 'LineStyle', ':')
xlabel('x')
ylabel('$\phi_{2}^{(1)}$')
legend show

figure()
plot(x1, real(phi12_B - phi12_con), 'DisplayName', 'r = 2')
title('difference between case B and control case')
xlabel('x')
ylabel('$\phi_{1}^{(2)}$')

figure()
plot(x2, real(phi23_B - phi23_con), 'DisplayName', 'r = 3')
hold on
title('difference between case B and control case')
plot(x2, real(phi23_B), 'DisplayName', 'Case B')
plot(x2, real(phi23_con), 'DisplayName', 'Case control')
plot(x2, flipud(real(phi23_B)), 'DisplayName', 'Case B flipped', 'LineStyle', '-.')
plot(x2, real( phi23_con - flipud(real(phi23_B))), 'DisplayName', 'Diff with Case B flipped', 'LineStyle', ':')
xlabel('x')
ylabel('$\phi_{2}^{(3)}$')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% CASE C: G1 and G2 coordinates swapped %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N1 = ceil(k*L1./(C_wl*2*pi)); N2 = ceil(k*L2./(C_wl*2*pi));% the number of boundary elements on screens 1 and 2

h1 = norm(r1endA-r1startA)/N1; % the length of each element on screen 1
h2 = norm(r2endB-r2startB)/N2; % the length of each element on screen 2
% Next calculate the x and y coordinates of the element midpoints on screen
% 1 and then on screen 2

x1 = r1startA(1)+((1:N1)-0.5)*(r1endA(1)-r1startA(1))/N1;
y1 = r1startA(2)+((1:N1)-0.5)*(r1endA(2)-r1startA(2))/N1;
x2 = r2startB(1)+((1:N2)-0.5)*(r2endB(1)-r2startB(1))/N2;
y2 = r2startB(2)+((1:N2)-0.5)*(r2endB(2)-r2startB(2))/N2;
h1vector = h1*ones(size(x1)); % Lengths of the elements on screen 1 
h2vector = h2*ones(size(x2)); % and on screen 2
h = [h1vector,h2vector]; % the lengths of the elements
x = [x1,x2]; % the x-coords of the element midpoints
y = [y1,y2]; % and their y-coords

phi1_C = bem2(x1,y1,h1vector,k,d);
phi2_C = zeros(size(x2));
phi_C = [phi1_C.',phi2_C];

phi10_C = phi1_C; %Keeping the zero iteration

phi2_C = bem2BS(x1,y1,x2,y2,h1vector,h2vector,phi1_C,k,d); % approximations to the values of phi at the element midppoints, phi(j) the value at (x(j),y(j))
phi_C = [phi1_C.',phi2_C.'];

%%% First iteration save
phi21_C = phi2_C; % keeping the 1 iteration


% Now let's compute phi12, the 2nd iteration
phi1_C = bem2BS(x2,y2,x1,y1,h2vector,h1vector,phi2_C,k,d);
phi_C = [phi1_C.',phi2_C.'];

phi12_C = phi1_C; 

% adapting from here on for more iterations

% Now let's compute phi23, the 3rd iteration
phi2_C = bem2BS(x1,y1,x2,y2,h1vector,h2vector,phi1_C,k,d); % approximations to the values of phi at the element midppoints, phi(j) the value at (x(j),y(j))
phi_C = [phi1_C.',phi2_C.'];

phi23_C = phi2_C;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%Comparison to control case%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
plot(x1, real(phi10_C - phi10_con), 'DisplayName', 'Difference')
hold on
plot(x1, real(phi10_C), 'DisplayName', 'Case C')
plot(x1, real(phi10_con), 'DisplayName', 'Case control')
plot(x1, flipud(real(phi10_C)), 'DisplayName', 'Case C flipped', 'LineStyle', '-.')
plot(x1, real( phi10_con - flipud(phi10_C) ), 'DisplayName', 'Difference, Case C flipped', 'LineStyle', ':')
title('difference between case C and control case')
xlabel('x')
ylabel('$\phi_{1}^{(0)}$')
legend show

figure()
plot(x2, real(phi21_C - phi21_con), 'DisplayName', 'Difference')
title('difference between case C and control case')
hold on
plot(x2, real(phi21_C), 'DisplayName', 'Case C')
plot(x2, real(phi21_con), 'DisplayName', 'Case control')
plot(x2, flipud(real(phi21_B)), 'DisplayName', 'Case C flipped', 'LineStyle', '-.')
plot(x2, real( phi21_con - flipud(real(phi21_B))), 'DisplayName', 'Diff with Case C flipped', 'LineStyle', ':')
xlabel('x')
ylabel('$\phi_{2}^{(1)}$')
legend show

figure()
plot(x1, real(phi12_C - phi12_con), 'DisplayName', 'Difference')
title('difference between case C and control case')
xlabel('x')
ylabel('$\phi_{1}^{(2)}$')
hold on
plot(x1, real(phi12_C), 'DisplayName', 'Case C')
plot(x1, real(phi12_con), 'DisplayName', 'Case control')
plot(x1, flipud(real(phi12_C)), 'DisplayName', 'Case C flipped', 'LineStyle', '-.')
plot(x1, real( phi12_con - flipud(phi12_C) ), 'DisplayName', 'Difference, Case C flipped', 'LineStyle', ':')
legend show

figure()
plot(x2, real(phi23_C - phi23_con), 'DisplayName', 'Difference')
title('difference between case C and control case')
hold on
plot(x2, real(phi23_C), 'DisplayName', 'Case C')
plot(x2, real(phi23_con), 'DisplayName', 'Case control')
plot(x2, flipud(real(phi23_C)), 'DisplayName', 'Case C flipped', 'LineStyle', '-.')
plot(x2, real( phi23_con - flipud(real(phi23_C))), 'DisplayName', 'Diff with Case C flipped', 'LineStyle', ':')
xlabel('x')
ylabel('$\phi_{2}^{(3)}$')
legend show




