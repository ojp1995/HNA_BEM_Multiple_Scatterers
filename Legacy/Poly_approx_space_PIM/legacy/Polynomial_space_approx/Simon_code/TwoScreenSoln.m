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
% and u_i is the incident plane wave
%
clear all
format compact
k = 5 % the wavenumber
r1start = [-2*pi, 2*pi]; r1end = [0, 0];
r2start = [4*pi, 0]; r2end = [7*pi, 3*pi];
L1 = sqrt( (r1end(2) - r1start(2))^2 + ( r1end(1) - r1start(1) )^2 )
L2 = sqrt( (r2end(2) - r2start(2))^2 + ( r2end(1) - r2start(1) )^2 )


C_wl = 1/640;

d = [0,-1]; % the direction of the incident wave
N1 = ceil(k*L1./(C_wl*2*pi)); N2 = ceil(k*L2./(C_wl*2*pi));% the number of boundary elements on screens 1 and 2

N = max(N1, N2);
N1 = N
N2 = N
h1 = norm(r1end-r1start)/N1; % the length of each element on screen 1
h2 = norm(r2end-r2start)/N2; % the length of each element on screen 2
% Next calculate the x and y coordinates of the element midpoints on screen
% 1 and then on screen 2
x1 = r1start(1)+((1:N1)-0.5)*(r1end(1)-r1start(1))/N1;
y1 = r1start(2)+((1:N1)-0.5)*(r1end(2)-r1start(2))/N1;
x2 = r2start(1)+((1:N2)-0.5)*(r2end(1)-r2start(1))/N2;
y2 = r2start(2)+((1:N2)-0.5)*(r2end(2)-r2start(2))/N2;
h = [h1*ones(size(x1)),h2*ones(size(x2))]; % the lengths of the elements
x = [x1,x2]; % the x-coords of the element midpoints
y = [y1,y2]; % and their y-coords
phi = bem2(x,y,h,k,d); % approximations to the values of phi at the element midppoints, phi(j) the value at (x(j),y(j))
phi1 = phi(1:N1); phi2 = phi(N1+1:end); % Split phi into the solution on screens 1 and 2
figure(1)
plot(((1:N1)-0.5)/N1,real(phi1))
xlabel('s/L_1')
ylabel('Re \phi_1(s)')
title('Real part of density on screen 1')
grid
figure(2)
plot(((1:N2)-0.5)/N2,real(phi2))
xlabel('s/L_2')
ylabel('Re \phi_2(s)')
title('Real part of density on screen 2')
grid

% Now plot the solution
X1 = min(x); X2 = max(x);
Y1 = min(y); Y2 = max(y);
dist = max(X2-X1,Y2-Y1);
Xlim(1) = X1-0.7*dist; Xlim(2) = X2+0.7*dist;
dist = 0.8*(Xlim(2)-Xlim(1))-(Y2-Y1);
Ylim(1) = round(Y1-dist/2); Ylim(2) = round(Y2 + dist/2);
lambda = 2*pi/k;
Np_perlamb = 10 % number of points per wavelength in the field plot in the domain
NX = ceil(Np_perlamb*(Xlim(2)-Xlim(1))/lambda)+1; % Mesh densities for the plots of the wave field
NY = ceil(Np_perlamb*(Ylim(2)-Ylim(1))/lambda)+1;
X = linspace(Xlim(1),Xlim(2),NX);
Y = linspace(Ylim(1),Ylim(2),NY);
[XX,YY] = meshgrid(X,Y);
u = ui(XX,YY,k,d(1),d(2)) + us(XX,YY,k,x,y,h,phi);
field_plot(3,XX,YY,u,'u',true)
hold on
plot([r1start(1),r1end(1)],[r1start(2),r1end(2)])
plot([r2start(1),r2end(1)],[r2start(2),r2end(2)])
hold off

SC_alin1_phi1_phi2_g1_n2pi_2pi_0_0_G2_4pi_0_7pi_3pi_Cwl640_d0_m1 = phi;

save('SC_alin1_phi1_phi2_g1_n2pi_2pi_0_0_G2_4pi_0_7pi_3pi_Cwl640_d0_m1')