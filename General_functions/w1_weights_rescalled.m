function w = w1_weights_rescalled(k, s, tj_1, tj, L)
% In this function we will be computing the integral which are the
% analytical weights for the PIM method for integrating the hankel function
% as a kernel.
%
%
% Inputs:
%
% Outputs:


w = -(s - L + tj_1).*(log( k*abs(s - L + tj_1) )) ...
    - (L - tj - s).*(log(k*abs( L - tj - s )))...
    + tj_1 - tj;

% now we need to find if there are any points where (s - tj_1) or (tj - s)
% is 0 and would cause a blow up and just override the weights with a - b

select1 = (s - tj_1) == 0;
select2 = (s - tj) == 0;

% if sum(select1 + select2)>0
%     keyboard
% end

w(select1) = -(2*tj_1(select1) - L).*log(k*(2*tj_1(select1) - L))...
    - (L - tj_1(select1) - tj(select1)).*log( k*(L - tj_1(select1) ...
    - tj(select1)) );
w(select2) = -(tj(select2) + tj_1(select2) - L).*log(k*(tj(select2) + ...
    tj_1(select2) - L)) + (L - 2*tj_1(select2)).*(log(k*(L - ...
    2*tj_1(select2))));



