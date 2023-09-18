function I = v1_weights_mid(a, b, k, x)
% In this function we will compute the weights needed for the PIM mid point
% rule when the singular part is ln|x-t|. We are computing the following
% integral for any x, I(x) = \int_{a}^{b} ln|x - t| dt.
%
%
% Problem parameters:
% a is the lower bound of the integral
% b is the upper bound of the integral
% x is a fixed point
%
% Assumptions:
% 1. a, b, x are scalar variable
% 2. Later will need to change for when x is far away from t.

% Due to the absolute values in the log we need to split approach the
% integrals in different ways:

% SEE BOTTOM OF CODE!!!
if a == x
    error('x cannot equal a yet')    
end

if b == x
    error('x cannot equal b yet')
end

if x < a  % when x < a, |x - t| = t - x
    
    I = (b - x)*log(b - x) + (x - a)*log(a - x) - (b - a);
    
elseif x > b  % when x > b, |x - t| = x - t
    
    I = (b - x)*log(x - b) + (x - a)*log(x - a) + a - b;
    
else  % when x \in [a, b].
   
    
    I = (x - a)*log(x - a) + (b - x)*log(b - x) + a - b;
    
end
%%%%
%%% We are actually integrating 
%%% \int_{a}^{b}ln(k|x - t|) dt =  \int_{a}^{b} (ln(k) + ln(|x - t|) ) dt
Ic = log(k)*(b - a); % this is ln(k) being separated out of the integral

I = I + Ic;  % combining both back together again.

end
    
    