function G_data = manipulate_collocation_points_graded(G_data)
% Manipulate collocation points so they are the arrange in the right
% orientation

G_data.s = [ G_data.t_mid_col(1:end) ; flip(G_data.L - ...
    G_data.t_mid_col(1:end)) ];
G_data.x_col = [ G_data.x_1_col(1:end) ; flip(G_data.x_2_col(1:end)) ];
G_data.y_col = [ G_data.y_1_col(1:end) ; flip(G_data.y_2_col(1:end)) ];
