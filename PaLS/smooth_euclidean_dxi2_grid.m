% Joshua Enxing
% Tufts University
% Done under supervision of Misha Kilmer and Eric Miller

% Finds the derivative of the smooth euclidean norm with respect to x2,
% basis function center

% Inputs:
%
% Meshgrid elements |X| and |Y|
% Parameter vector |p|
% Nu value |v| that is used with the smooth euclidean norm

% Outputs:
%
% Matrix of derivative of smooth euclidean norm with respect to x2 on the grid
function N = smooth_euclidean_dxi2_grid(X,Y,p,v)

j = length(p)/4;
N = zeros(length(X));
for i = 1:j
    A = p(4*(i-1)+3)*ones(length(X));
    B = p(4*(i-1)+4)*ones(length(Y));
    C = p(4*(i-1)+2)*(X - A);
    D = p(4*(i-1)+2)*(Y - B);
    R = (C.^2+D.^2+v^2*ones(length(X))).^0.5;
    N = N + (p(4*(i-1)+2)^2*(p(4*(i-1)+4)*ones(length(X)) - Y))./R;
end