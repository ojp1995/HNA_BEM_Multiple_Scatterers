function us = HF_it_soln_in_2D(G1_data, phi1, G2_data, phi2, k, X, Y)
% This function will compute the scattered field in 2D

y1q1 = G1_data.x_q_comb_outer;
y1q2 = G1_data.y_q_comb_outer;
w1 = G1_data.w_comb_outer;

y2q1 = G2_data.x_q_comb_outer;
y2q2 = G2_data.y_q_comb_outer;
w2 = G2_data.w_comb_outer;

us = zeros(length(Y), length(X));
for ix  = 1:length(X)

    for iy = 1:length(Y)

        dist1 = sqrt( (X(ix) - y1q1).^2 + (Y(iy) - y1q2).^2 );

        dist2 = sqrt( (X(ix) - y2q1).^2 + (Y(iy) - y2q2).^2 );

        us(iy, ix) = 1i*sum(phi1.*w1.*besselh(0, k*dist1))/4 ...
            + 1i*sum(phi2.*w2.*besselh(0, k*dist2))/4;

    end 

end
