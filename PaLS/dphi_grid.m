% Joshua Enxing
% Tufts University
% Done under supervision of Misha Kilmer and Eric Miller

% Evaluates derivative of phi on the grid with input R

% Inputs:
%
% Matrix of values |R|

% Outputs:
% 
% Matrix |dphi| of derivative of phi on the grid
function dphi = dphi_grid(R)

dphi = R < 1;
dphi = R.*dphi;
dphi = -20.*dphi.*(1-dphi).*(1-dphi).*(1-dphi);