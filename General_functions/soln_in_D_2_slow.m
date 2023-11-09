function us = soln_in_D_2_slow(G_data1, phi1, G_data2, phi2, k, X, Y)
% Computing the solution in the domain slowly using a double for loop,
% takes much longer
%
% G_data, all the information about the screen
% phi, solution evaluated at quadrature nodes
% k, wavenumber
% X, Xmesh - vector
% Y, Ymesh - vector

y1q1 = [G_data1.x_1_q ; flip(G_data1.x_2_q)];
y1q2 = [G_data1.y_1_q ; flip(G_data1.y_2_q)];
w1 = [G_data1.w ; flip(G_data1.w)];

y2q1 = [G_data2.x_1_q ; flip(G_data2.x_2_q)];
y2q2 = [G_data2.y_1_q ; flip(G_data2.y_2_q)];
w2 = [G_data2.w ; flip(G_data2.w)];

us = zeros(length(Y), length(X));


for ix = 1:length(X)

    for iy = 1:length(Y)

        dist1 = sqrt( (X(ix) - y1q1).^2 + (Y(iy) - y1q2).^2 );

        dist2 = sqrt( (X(ix) - y2q1).^2 + (Y(iy) - y2q2).^2 );

        us(iy, ix) = 1i*sum(phi1.*w1.*besselh(0, k*dist1))/4 ...
            + 1i*sum(phi2.*w2.*besselh(0, k*dist2))/4;

    end
end
