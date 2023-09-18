function [S11, S12, S21, S22, u_inc1, uinc2] = ...
    compute_matrices_for_iterative_solve(G1, G2, k, )
% In this function we will be computing the 4 matrices S11, S12, S21 and
% S22 to either be used for the all in one or iterative solver, the
% corresponding right hand side vector will also be given for the initial
% incident case.
%
% The matrices will be computed with piecewise constant basis functions
% where the intervals they are supported on will be given as an input.
%
% Input parameters:
% G1 - Coordinates for Gamma1
% G2 - coordinates for Gamma2
% k wavenumber


% C1, C2, constants for product integration method

% t1_grid - Grid of support for basis functions for Gamma1
% t1_grid - Grid of support for basis functions for Gamma2


for n = 1:length(t1_grid)  % basis function loop
    % computation of S11
    for m = 1:length(x1_col)  % collocation loop

        S11(m, n) = 1;  % holding for structure

    end

    % computation of S21
    for m = 1:length(x2_col) % collocation loop for S21

        S21(m, n) = 0;  % holding for structure

    end 

end

for n = 1:length(t2_grid)  % basis function loop
    % computation of S22
    for m = 1:length(x2_col)  % collocation loop

        S22(m, n) = 1;  % holding for structure

    end

    % computation of S12
    for m = 1:length(x1_col) % collocation loop for S12

        S12(m, n) = 0;  % holding for structure

    end 

end

% computing incident for each

u_inc1 = 'holding';
u_inc2 = 'holding'




