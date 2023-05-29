%
% Computation by a piecewise constant collocation method of the solution 
% phi, on a screen Gamma consisting of two straight lines, by solving
%
% S phi = u^i|_Gamma,
% 
% where S is the single-layer potential integral operator, with kernel
%
% Phi(r,r_0) = (i/4) H_0^{(1)}(k|r-r_0|)
%
% and u_i is the incident plane wave.
%
% This file computes the first two iterations of our multiple scattering
% solution and the corresponding field.
%
format compact
k = 5 % the wavenumber
% case 1 \ /
% r1start = [-2*pi,2*pi]; r1end = [0,0]; % Endpoints of screen 1
% r2start = [4*pi,0]; r2end = [6*pi,2*pi]; %Endpoints of screen 2

%case 2 \ / but bigger
r1start = [-2*pi, 2*pi]; r1end = [0, 0];
r2start = [4*pi, 0]; r2end = [7*pi, 3*pi];
L1 = sqrt( (r1end(2) - r1start(2))^2 + ( r1end(1) - r1start(1) )^2 )
L2 = sqrt( (r2end(2) - r2start(2))^2 + ( r2end(1) - r2start(1) )^2 )

%case 3 \|, closer together, more reflections
 C_wl = 1/40;

d = [0,-1]; % the direction of the incident wave
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
% The next line does the zero iteration, approximating the values of phi1 
% at the element midpoints, phi(j) the value at (x(j),y(j))
phi1 = bem2(x1,y1,h1vector,k,d);
phi2 = zeros(size(x2));
phi = [phi1.',phi2];
figure()
plot(((1:N1)-0.5)/N1,real(phi1))
xlabel('s/L_1')
ylabel('Re \phi_1^{(0)}(s)')
title('Real part of density on screen 1 after 0 iteration')
grid
ut_r0 = ComputeAndPlotTotalField(x,y,h,k,d,phi,'u after iteration 0',...
    r1start,r1end,r2start,r2end);
disp('Pausing after iteration zero. Press any key to continue ...')

%%%% zeroth iteration save
phi10 = phi1; %Keeping the zero iteration

phi2 = bem2BS(x1,y1,x2,y2,h1vector,h2vector,phi1,k,d); % approximations to the values of phi at the element midppoints, phi(j) the value at (x(j),y(j))
phi = [phi1.',phi2.'];
figure()
plot(((1:N2)-0.5)/N2,real(phi2))
xlabel('s/L_2')
ylabel('Re \phi_2^{(1)}(s)')
title('Real part of density on screen 2 after iteration 1')
grid
ut_r1 = ComputeAndPlotTotalField(x,y,h,k,d,phi,'u after iteration 1',...
    r1start,r1end,r2start,r2end);
disp('Pausing after iteration one. Press any key to continue ...')

%%% First iteration save
phi21 = phi2; % keeping the 1 iteration


% Now let's compute phi12, the 2nd iteration
phi1 = bem2BS(x2,y2,x1,y1,h2vector,h1vector,phi2,k,d);
phi = [phi1.',phi2.'];
figure()
plot(((1:N1)-0.5)/N1,real(phi10),((1:N1)-0.5)/N1,real(phi1))
xlabel('s/L_1')
ylabel('Re \phi_1^{(r)}(s)')
legend('r = 0','r = 2')
title('Real part of density on screen 1')
grid
ut_r2 = ComputeAndPlotTotalField(x,y,h,k,d,phi,'u after iteration 2',...
    r1start,r1end,r2start,r2end);

%%% 2nd iteration save
phi12 = phi1; 

% adapting from here on for more iterations

% Now let's compute phi23, the 3rd iteration
phi2 = bem2BS(x1,y1,x2,y2,h1vector,h2vector,phi1,k,d); % approximations to the values of phi at the element midppoints, phi(j) the value at (x(j),y(j))
phi = [phi1.',phi2.'];
figure()
plot(((1:N2)-0.5)/N2,real(phi21),((1:N2)-0.5)/N2,real(phi2))
xlabel('s/L_1')
ylabel('Re \phi_2^{(r)}(s)')
legend('r = 1','r = 3')
title('Real part of density on screen 2')
grid
ut_r3 = ComputeAndPlotTotalField(x,y,h,k,d,phi,'u after iteration 3',...
    r1start,r1end,r2start,r2end, 4);

% 3rd iteration save
phi23 = phi2;

% saving data:
% SCphi10_g1_n2pi_2pi_0_0_G2_4pi_0_7pi_3pi_Cwl640_d0_m1 = phi10;
% save('SCphi10_g1_n2pi_2pi_0_0_G2_4pi_0_7pi_3pi_Cwl640_d0_m1')
% 
% SCphi21_g1_n2pi_2pi_0_0_G2_4pi_0_7pi_3pi_Cwl640_d0_m1 = phi21;
% save('SCphi21_g1_n2pi_2pi_0_0_G2_4pi_0_7pi_3pi_Cwl640_d0_m1')
% 
% SCphi12_g1_n2pi_2pi_0_0_G2_4pi_0_7pi_3pi_Cwl640_d0_m1 = phi12;
% save('SCphi12_g1_n2pi_2pi_0_0_G2_4pi_0_7pi_3pi_Cwl640_d0_m1')
% forgot I only need to save one lot of data and the rest follows.

newSCphi23_g1_n2pi_2pi_0_0_G2_4pi_0_7pi_3pi_Cwl640_d0_m1 = phi23;
save('newSCphi23_g1_n2pi_2pi_0_0_G2_4pi_0_7pi_3pi_Cwl640_d0_m1')
