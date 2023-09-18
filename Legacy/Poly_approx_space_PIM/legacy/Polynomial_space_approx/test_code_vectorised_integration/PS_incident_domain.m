function z = PS_incident_domain(xps, yps, X, Y, k)
% computing the point source in the domian.

% (xps, yps) is the coordinate of the point source
% X is the vector in the domain in the X direction
% Y is the vector in the domain in the y direction
% k is the wavenumber 

for yj = 1:length(Y)
    for xj = 1:length(X)
        % prevent blow up at point source
               if (X(xj) ~= xps) && ( Y(yj) ~= yps)
                   z(yj, xj) = besselh(0, k*sqrt( (X(xj) - xps).^2 + ( Y(yj) - yps ).^2 ));
               else
                   z(yj, xj) = 0;
               end
       
    end
end