% Joshua Enxing
% Tufts University
% Done under supervision of Misha Kilmer and Eric Miller

% Finds the value of the derivative of the smooth heaviside function on the grid

% Inputs:
%
% Matrix of values want to evaluate at, |M|
% Epsilon value |eps| that is used with the smooth heaviside function

% Outputs:
%
% Matrix of derivative of smooth heaviside values |V| on the grid
function V = smooth_heaviside_deriv_grid(M,eps)

L = abs(M) <= eps;
V = M.*L;
V = (1/(2*eps))*(1+cos((pi*(V))/eps));
V = V.*L;
