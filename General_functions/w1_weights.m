function w = w1_weights(k, s, tj_1, tj)
% In this function we will be computing the integral which are the
% analytical weights for the PIM method for integrating the hankel function
% as a kernel.
%
%
% Inputs:
%
% Outputs:


% problem in this function that is causing blow up, sloppy way to get
% around it:
% if sum((tj - s)==0) || sum(((s - tj_1))==0)
%     for j = 1:length(tj_1)
%         if s == tj_1(j) || s==t_j(j)
%             w(j, 1) = 0


w = (s - tj_1).*log( k*abs(s - tj_1) ) + (tj - s).*log(k*abs( tj - s )) + tj_1 - tj;

% now we need to find if there are any points where (s - tj_1) or (tj - s)
% is 0 and would cause a blow up and just override the weights with 0.
% 
% lower_zero = ( (s - tj_1) == 0);
% upper_zero = (   (s - tj) == 0);
% 

if (sum((s - tj_1) == 0))> 0 || (sum((s - tj) == 0)>0)
    lower_index = find((s - tj_1) == 0, 10);
    upper_index = find((s - tj) == 0, 10);
    
    w(lower_index) = tj_1(lower_index) - tj(lower_index);
    
    w(upper_index) = tj_1(upper_index) - tj(upper_index);

end

% 
% 
% w(lower_index) = 0;
% 
% w(upper_index) = 0;


