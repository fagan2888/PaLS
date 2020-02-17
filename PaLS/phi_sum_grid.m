% Joshua Enxing
% Tufts University
% Done under supervision of Misha Kilmer and Eric Miller

% Evaluates phi, sum of basis functions, on the given grid

% Inputs:
% 
% Parameter vector |p|
% Meshgrid elements |X| and |Y|
% Nu value |v| that is used with the smooth euclidean norm
% Option |alpha_yes| for if alpha is used or not, alpha_yes=0 means alpha
% is not used

% Outputs:
% 
% Matrix of phi values on entire grid |phi|, matrix |R| that is recycled
% for computational efficiency
function [phi,R] = phi_sum_grid(p,X,Y,v,alph_yes)
  
%Number of basis functions j
j = length(p)/4;
R = zeros(length(X));
phi = R;

%a is vector of amplitudes
a = p(1:4:end);

%b is vector of dilation factors
b = p(2:4:end);

%xx is vector of x centers
xx = p(3:4:end);

%yy is vector of y centers
xy = p(4:4:end);

%Iterate over all basis functions
for i=1:j
    
    %Vectorization of equation
    L = b(i)*(X - xx(i));
    D = b(i)*(Y - xy(i));
    R = sqrt(L.*L+D.*D+v^2);
    L = max(0,1-R);
    
    %For computational efficiency 
    L = L.*L;
    L = L.*L;
    if alph_yes == 0
        phi = phi + (L.*(4*R+1));
    else
        phi = phi + a(i)*(L.*(4*R+1));
    end
end