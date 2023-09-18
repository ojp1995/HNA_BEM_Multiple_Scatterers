function u = ui(x,y,k,dx,dy)
%
% function u = ui(x,y,k,dx,dy)
%
% This function calculates a time harmonic plane wave solution of the
% Helmholtz equation.
%
% The inputs are:
%
% x,y     real
%         (x,y) are the coordinates of the position vector at which to
%         compute the plane wave.
%
% k       real or complex
%         is the wave number, with Re k > 0, Im k >= 0.
%
% dx, dy  real
%         (dx,dy) are the components of the unit vector in the direction of
%         the plane wave.
%
%
% The output is:
% 
% u       complex
%         The value of the plane wave at (x,y).
%
% N.B. If x,y are vectors of the same length then u will also be a vector,
% with u(j) the value of the plane wave at the point (x(j),y(j)).
% Similarly, if x and y are matrices of the same size then u will be a
% matrix of the same size.
%
u = exp(i*k*(x*dx+y*dy));