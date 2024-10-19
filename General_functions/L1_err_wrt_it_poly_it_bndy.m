function err = L1_err_wrt_it_poly_it_bndy(G_data, err_quad, aj_R, R_max)
% In this function we will compute the L1 error for solutions of the
% iteratirve method with respect to the NUMBER OF ITERATIONS. The true
% solution will be taken as the largest number of iterations.
% 
% Comparing solution on BOUNDARY
%
% Inputs:
% G_data, all the information about the G1, most important is the bf
% discretisation
% err_quad, object containing all info about the err_quad
% aj_R, coefficients for Gamma.
% R_max, maximum number of iterations

% Outputs:
% err, vector of normalised L1 error

% first computing the solution on the bndy
phi1 = zeros(R_max, 2*length(err_quad.x_1_q));
for r = 1:R_max
    
    phi1(r, :) = graded_coeff_2_solution(aj_R(:, r), G_data.t_bf_grid, ...
        err_quad.t_mid_q, G_data.L);

end

% Now computing the L1 error
err = zeros(R_max - 1, 1);
for r = 1:R_max - 1

%     err(r) = (abs(phi1(end, :) - phi1(r, :))./abs(phi1(end, :)))...
%         *[err_quad.w ; flip(err_quad.w)];

    err(r) = (abs(phi1(end, :) - phi1(r, :))*[err_quad.w ; flip(err_quad.w)])...
        ./(abs(phi1(end, :)))*[err_quad.w ; flip(err_quad.w)];
end


