function G2_data = get_coeff_for_different_sides_of_screen(G1_data, ...
    G2_data)
% In this function we will determine the parameter beta, which decides
% whether G2 is incident on the positive side or negative side of Gamma1

% We are going to compute the intersections of the screens and depending on
% that we will give a negative, positive or mixture depending on how the
% screens intersect

% first computing the gradient
try

    m1 = (G1_data.G(4) - G1_data.G(2))/(G1_data.G(3) - G1_data.G(1));
    m2 = (G2_data.G(4) - G2_data.G(2))/(G2_data.G(3) - G2_data.G(1));

catch
    G1_x_diff = G1_data.G(3) - G1_data.G(1)
    G2_x_diff = G2_data.G(3) - G2_data.G(1)
    disp('Attempted divide by 0?')

end

% Extension of Gamma_{1} extending into Gamma_{2} or the extension of it

x_intersect = (-m2*G2_data.G(3) + G2_data.G(4) + m1*G1_data.G(3) ...
    - G1_data.G(4))/(m1 - m2);

if G2_data.G(1) < x_intersect && x_intersect < G2_data.G(3)

    positive_side_outer = zeros(length(G2_data.x_q_comb_outer), 1);
    negative_side_outer = zeros(length(G2_data.x_q_comb_outer), 1);

    positive_side_inner = zeros(length(G2_data.x_q_comb_inner), 1);
    negative_side_inner = zeros(length(G2_data.x_q_comb_inner), 1);

    positive_side_outer = (x_intersect <= G2_data.x_q_comb_outer); 
    negative_side_outer = (x_intersect > G2_data.x_q_comb_outer);

    G2_data.beta_outer = positive_side_outer - negative_side_outer;

    positive_side_inner = (x_intersect <= G2_data.x_q_comb_inner); 
    negative_side_inner = (x_intersect > G2_data.x_q_comb_inner);

    G2_data.beta_inner = positive_side_inner - negative_side_inner;

else  %computing which side is which

    if min(G2_data.G(1), G2_data.G(3)) > max(G1_data.G(1), G1_data.G(3)) % Gamma2 is to the right of Gamma1
        
        if x_intersect < max(G2_data.G(1), G2_data.G(3))
        % negative case, Gamma2 is on the negative side of Gamma1
        G2_data.beta_outer = ones(length(G2_data.x_q_comb_outer), 1);
        G2_data.beta_inner = ones(length(G2_data.x_q_comb_inner), 1);

        else
            % positive side
            G2_data.beta_outer = -ones(length(G2_data.x_q_comb_outer), 1);
            G2_data.beta_inner = -ones(length(G2_data.x_q_comb_inner), 1);
    
        end

    else % Gamma2 is to the left of Gamma1


        if x_intersect < max(G2_data.G(1), G2_data.G(3))
            % negative case, Gamma2 is on the negative side of Gamma1
            G2_data.beta_outer = -ones(length(G2_data.x_q_comb_outer), 1);
            G2_data.beta_inner = -ones(length(G2_data.x_q_comb_inner), 1);
    
        else
            % positive side
            G2_data.beta_outer = ones(length(G2_data.x_q_comb_outer), 1);
            G2_data.beta_inner = ones(length(G2_data.x_q_comb_inner), 1);
    
        end

    end

end


% sanity check

if (sum(abs(G2_data.beta_outer)) ~= length(G2_data.x_q_comb_outer))
    error('Length for outer beta coeff is different')
end

if (sum(abs(G2_data.beta_inner)) ~= length(G2_data.x_q_comb_inner))
    error('Length for inner is different')
end

