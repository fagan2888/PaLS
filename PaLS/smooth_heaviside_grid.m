% Joshua Enxing
% Tufts University
% Done under supervision of Misha Kilmer and Eric Miller

% Finds the value of the smooth heaviside function on the grid

% Inputs:
%
% Matrix of values want to evaluate at, |M|
% Epsilon value |eps| that is used with the smooth heaviside function

% Outputs:
%
% Matrix of smooth heaviside values |V| on the grid
function V = smooth_heaviside_grid(M,eps)

V = zeros(length(M));
V(M > eps) = 1;
L = abs(M) <= eps;
A = 0.5*(ones(length(M))+(M./eps)+(1/pi).*sin((pi.*M)./eps));
B = A.*L;
V = V + B;
