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

elseif m1 == m2  % parallel screen case

    if G2_data.G(1) == G2_data.G(3)  % parallel x (similar to | |)

        if G2_data.G(2) > G1_data.G(2)  % Gamma2 above Gamma1
            G2_data.beta_outer = ones(length(G2_data.x_q_comb_outer), 1);
            G2_data.beta_inner = ones(length(G2_data.x_q_comb_inner), 1);
        else
            % Gamma2 below Gamma1
            G2_data.beta_outer = -ones(length(G2_data.x_q_comb_outer), 1);
            G2_data.beta_inner = -ones(length(G2_data.x_q_comb_inner), 1);
        end

    elseif G2_data.G(2) == G2_data.G(4)  % parallel y (similar to =)
        if G2_data.G(2) > G1_data.G(2)  % Gamma2 to the right of (facing front) Gamma1
            G2_data.beta_outer = ones(length(G2_data.x_q_comb_outer), 1);
            G2_data.beta_inner = ones(length(G2_data.x_q_comb_inner), 1);
        else
            % Gamma2 to the left of (behind) Gamma1
            G2_data.beta_outer = -ones(length(G2_data.x_q_comb_outer), 1);
            G2_data.beta_inner = -ones(length(G2_data.x_q_comb_inner), 1);
        end

    else  % case that parallel off angle
        error('Parallel case, not aligned to x or y axis, not computed yet.')

    end

else  % in the case that it is either all acting on positive or negative
    
    % Both screens are to the left of the intersection
    if min(G2_data.G(1), G2_data.G(3)) < x_intersect &&  ...
            min(G1_data.G(1), G1_data.G(3)) < x_intersect 

        % now we want to determine is Gamma2 facing the positive or
        % negative side of Gamma1

        if min(G2_data.G(2), G2_data.G(4)) < max(G1_data.G(2), G1_data.G(4)) % Gamma 2 below Gamma1

            G2_data.beta_outer = -ones(length(G2_data.x_q_comb_outer), 1);
            G2_data.beta_inner = -ones(length(G2_data.x_q_comb_inner), 1); 

        else
            % Gamma2 is above Gamma1
            G2_data.beta_outer = ones(length(G2_data.x_q_comb_outer), 1);
            G2_data.beta_inner = ones(length(G2_data.x_q_comb_inner), 1);
        end

        


    % Both screens are to the right
    elseif min(G2_data.G(1), G2_data.G(3)) > x_intersect &&  ...
            min(G1_data.G(1), G1_data.G(3)) > x_intersect 
        
        % now we want to determine is Gamma2 facing the positive or
        % negative side of Gamma1

        % is Gamma 2 below Gamma1
        if min(G2_data.G(2), G2_data.G(4)) < max(G1_data.G(2), G1_data.G(4))
            
            G2_data.beta_outer = -ones(length(G2_data.x_q_comb_outer), 1);
            G2_data.beta_inner = -ones(length(G2_data.x_q_comb_inner), 1); 

        else
            % Gamma2 is above Gamma1
            G2_data.beta_outer = ones(length(G2_data.x_q_comb_outer), 1);
            G2_data.beta_inner = ones(length(G2_data.x_q_comb_inner), 1);

        end


    % Screens are facing each other
    else

        if dot(G2_data.n, G1_data.n) <= 0
            % screens are facing each other
            G2_data.beta_outer = ones(length(G2_data.x_q_comb_outer), 1);
            G2_data.beta_inner = ones(length(G2_data.x_q_comb_inner), 1);

        else
            % screens are not facing each other
            warning('Not sure how this cases exists!')

        end

    end

end

%     if min(G2_data.G(1), G2_data.G(3)) > max(G1_data.G(1), G1_data.G(3)) % Gamma2 is to the right of Gamma1
%         
%         if x_intersect < max(G2_data.G(1), G2_data.G(3))
%         % negative case, Gamma2 is on the negative side of Gamma1
%         G2_data.beta_outer = ones(length(G2_data.x_q_comb_outer), 1);
%         G2_data.beta_inner = ones(length(G2_data.x_q_comb_inner), 1);
% 
%         else
%             % positive side
%             G2_data.beta_outer = -ones(length(G2_data.x_q_comb_outer), 1);
%             G2_data.beta_inner = -ones(length(G2_data.x_q_comb_inner), 1);
%     
%         end
% 
%     else % Gamma2 is to the left of Gamma1
% 
% 
%         if x_intersect < max(G2_data.G(1), G2_data.G(3))
%             % negative case, Gamma2 is on the negative side of Gamma1
%             G2_data.beta_outer = -ones(length(G2_data.x_q_comb_outer), 1);
%             G2_data.beta_inner = -ones(length(G2_data.x_q_comb_inner), 1);
%     
%         else
%             % positive side
%             G2_data.beta_outer = ones(length(G2_data.x_q_comb_outer), 1);
%             G2_data.beta_inner = ones(length(G2_data.x_q_comb_inner), 1);
%     
%         end
% 
%     end
% 
% end
% 
% 
% % sanity check
% 
% if (sum(abs(G2_data.beta_outer)) ~= length(G2_data.x_q_comb_outer))
%     error('Length for outer beta coeff is different')
% end
% 
% if (sum(abs(G2_data.beta_inner)) ~= length(G2_data.x_q_comb_inner))
%     error('Length for inner is different')
% end
% 
