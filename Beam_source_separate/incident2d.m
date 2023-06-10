function ui = incident2d(k, theta, X, Y)
% Compute the incident field in 2D

ui = zeros(length(Y), length(X));

for ix = 1:length(X)
    for iy = 1:length(Y)
        
        ui(iy, ix) = incident(k, theta, X(ix), Y(iy));
        
    end 
end