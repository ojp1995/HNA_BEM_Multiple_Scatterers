function ui = incident_domain(k, theta, X, Y)
% This function will compute the incident field in the domain. 
% ui = exp(ikx\cdot d), x is the point in the domain and d is the direction
% of the wave d = (sin(theta), -cos(theta)), theta is the angle between the
% downwards vertical and the incident wave anti clockwise
%
%
% Problem parameters:
% k is the wavenumber
% theta is the angle between the dowards vertical and the incident field
%
% Discretisation parameters
% X is the discretisation in the X direction
% Y is the discretisation in the Y direction

ui = zeros(length(Y), length(X));

for ix = 1:length(X)% loop through every point in x direction
    
    for iy = 1:length(Y)  % looops through every point in the y direction
        
        ui(iy, ix) = exp(1i*k*(X(ix)*sin(theta) - Y(iy)*cos(theta)));
        
    end
end