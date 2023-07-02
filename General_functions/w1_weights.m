function w = w1_weights(k, s, tj_1, tj)
% In this function we will be computing the integral which are the
% analytical weights for the PIM method for integrating the hankel function
% as a kernel.
%
%
% Inputs:
%
% Outputs:


w = (s - tj_1).*(log( k*abs(s - tj_1) )) + (tj - s).*(log(k*abs( tj - s )))...
    + tj_1 - tj;

% now we need to find if there are any points where (s - tj_1) or (tj - s)
% is 0 and would cause a blow up and just override the weights with a - b

select1 = (s - tj_1) == 0;
select2 = (s - tj) == 0;

% if sum(select1 + select2)>0
%     keyboard
% end

w(select1) = (tj(select1) - tj_1(select1)).*log(k*(tj(select1) ...
    - tj_1(select1))) + tj_1(select1) - tj(select1);
w(select2) = (tj(select2) - tj_1(select2)).*log(k*(tj(select2) ...
    - tj_1(select2))) + tj_1(select2) - tj(select2);



