function err = HF_it_L1_err_wrt_R(phi_true, phi_test, R_max, err_quad)
% computing L1 error for HF_it method

err = zeros(R_max-1, 1);

for r = 1:R_max - 1

%     err(r) = sum((abs(phi_true - phi_test{r})./abs(phi_true)).*err_quad);
    err(r) = sum((abs(phi_true - phi_test{r}).*err_quad))/ ...
        sum(abs(phi_true).*err_quad);

end