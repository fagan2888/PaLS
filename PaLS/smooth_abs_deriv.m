% Joshua Enxing
% Tufts University
% Done under supervision of Misha Kilmer and Eric Miller

% Derivative of smooth absolute value function

% Inputs:
%
% Object you want to find derivative of |s|
% Epsilon value |eps|

% Ouputs:
%
% |new|, smooth abs derivative value of |s| 
function der  = smooth_abs_deriv(s,eps)

new = smooth_abs(s,eps);
der = s./new;