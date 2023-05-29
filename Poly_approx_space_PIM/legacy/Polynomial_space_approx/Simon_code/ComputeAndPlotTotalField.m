function u = ComputeAndPlotTotalField(x,y,h,k,d,phi,Title,r1start,r1end,...
                r2start,r2end,FigNo,Xlim,Ylim)
% function ComputeAndPlotTotalField(x,y,h,k,d,phi,Title,r1start,r1end,...
%                r2start,r2end,FigNo,Xlim,Ylim)
%
% This function plots the total field in the domain, assuming that
% is represented as the sum of an incident plane wave and a scattered
% field that is the field due to a single layer potential (SLP) with 
% density phi supported on two screens. The method is to appproximate the 
% SLP by a simple composite midpoint rule.
%
% The inputs are:
%
% x,y,h      real vectors of the same length N
%            The screens on which the SLP is supported are divided into 
%            N elements in total. (x(i),y(i)) are the 
%            coordinates of the midpoint and h(i) the length of the ith
%            element, i=1,2,...,N
%
% k          real positive
%            The wave number
%
% d          real vector of length 2
%            A unit vector in the direction of the plane wave
%
% phi        complex vector of length N
%            phi(i) is the density of the SLP at (x(i),y(i)), i=1,...,N.
%
% Title      string
%            The title for the plot
%
% r1start    real vector length 2
%            The coordinates of one end of screen 1
%
% r1end      real vector length 2
%            The coordinates of the other end of screen 1
%
% r2start    real vector length 2
%            The coordinates of one end of screen 2
%
% r2end      real vector length 2
%            The coordinates of the other end of screen 2
%
% FigNo      integer, positive
%            The figure number to be used
%            Default is 3
%
% Xlim,Ylim  real vectors of length 2
%            The field is plotted for x between Xlim(1) and Xlim(2), for
%            y between Ylim(1) and Ylim(2).
%            If these arguments are missing, the function computes sensible
%            defaults.
%
if nargin < 12
    FigNo = 3;
end
if nargin < 13    
    X1 = min(x); X2 = max(x);
    Y1 = min(y); Y2 = max(y);
    dist = max(X2-X1,Y2-Y1);
    Xlim(1) = X1-0.7*dist; Xlim(2) = X2+0.7*dist;
    dist = 0.8*(Xlim(2)-Xlim(1))-(Y2-Y1);
    Ylim(1) = round(Y1-dist/2); Ylim(2) = round(Y2 + dist/2);
end
lambda = 2*pi/k;
Np_perlamb = 10; % number of points per wavelength in the field plot in the domain
NX = ceil(Np_perlamb*(Xlim(2)-Xlim(1))/lambda)+1; % Mesh densities for the plots of the wave field
NY = ceil(Np_perlamb*(Ylim(2)-Ylim(1))/lambda)+1;
X = linspace(Xlim(1),Xlim(2),NX);
Y = linspace(Ylim(1),Ylim(2),NY);
[XX,YY] = meshgrid(X,Y);
u = ui(XX,YY,k,d(1),d(2)) + us(XX,YY,k,x,y,h,phi);
field_plot(FigNo,XX,YY,u,Title,true)
hold on
plot([r1start(1),r1end(1)],[r1start(2),r1end(2)])
plot([r2start(1),r2end(1)],[r2start(2),r2end(2)])
hold off
end

