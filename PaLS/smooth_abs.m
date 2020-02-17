% Joshua Enxing
% Tufts University
% Done under supervision of Misha Kilmer and Eric Miller

% Smooth version of absolute value function

% Inputs:
%
% Object you want to find abs value of |s|
% Epsilon value |eps|

% Ouputs:
%
% |new|, smooth abs value of |s| 

function new = smooth_abs(s,eps)

new = zeros(size(s));
for i=1:length(s)
    new(i) = sqrt(s(i)^2+eps);
end