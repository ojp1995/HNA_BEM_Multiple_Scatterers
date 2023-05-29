function S21 = S21_op_robhop(x1, y1, x2, y2, k, h)
% REWRITTEN as robust hopefully.
% This function computes S21\phi_{1}^{(r)} = \int_{\G1}\phik(x,
% y)\phi_{1}^{(r)}(y) ds(y), x \in \G2
% Idea here is the beam originating from G1 incident on G2.
%
% Can produce S12 op if coordinates are flipped at input.
%
% Problem parameters
% (x1, y1) is the coordinates of G1
% (x2, y2) are the coordinates of G2
% k is the wavenumber
%
% Discretisation parameters:
% h are the integration weights

for j = 1: length(x1)  %looping through the integration intervals on \Gamma_1
    
    for l = 1:length(x2)  % looping through the collocation points on \Gamma_2
        dist = sqrt( (x2(l) - x1(j))^2 + (y2(l) - y1(j))^2 );
        S21(l, j) = 1i*h(j)*besselh(0, k*dist)/4;
        
    end
    
end