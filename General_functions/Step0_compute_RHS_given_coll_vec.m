function f1_0 = Step0_compute_RHS_given_coll_vec(kwave, theta, ...
    col_points1, x1_col, y1_col, L1, vertices1, d, h, nq, tq, C1, C2 )

u_inc = incident(kwave, theta, x1_col, y1_col);

f_duidn = duidn(vertices1, L1, kwave, d, nq);
LoB = 2*PIM_int_hankel_f(kwave, col_points1, h, nq, f_duidn, tq, C1, C2);

f1_0 = u_inc.' - LoB;