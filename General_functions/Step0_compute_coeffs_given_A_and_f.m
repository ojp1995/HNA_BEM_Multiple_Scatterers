function v_N_G1_0 = Step0_compute_coeffs_given_A_and_f(A1, f1_0, VHNA1)

aj = A1\f1_0; % computes the coefficients

v_N_G1_0 = ProjectionFunction(aj, VHNA1);