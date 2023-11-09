function us = soln_in_D_slow(G_data, phi, k, X, Y)
% Computing the solution in the domain slowly using a double for loop,
% takes much longer
%
% G_data, all the information about the screen
% phi, solution evaluated at quadrature nodes
% k, wavenumber
% X, Xmesh - vector
% Y, Ymesh - vector

yq1 = [G_data.x_1_q ; flip(G_data.x_2_q)];
yq2 = [G_data.y_1_q ; flip(G_data.y_2_q)];
w = [G_data.w ; flip(G_data.w)];

us = zeros(length(Y), length(X));
% us = zeros(size(X));

for ix = 1:length(X)

    for iy = 1:length(Y)

        dist = sqrt( (X(ix) - yq1).^2 + (Y(iy) - yq2).^2 );

        us(iy, ix) = 1i*sum(phi.*w.*besselh(0, k*dist))/4;

    end
end
