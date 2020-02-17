% Joshua Enxing
% Tufts University
% Done under supervision of Misha Kilmer and Eric Miller

% Evaluates f on the given grid

% Inputs:
% 
% Parameter vector |p|
% Meshgrid elements |X| and |Y|
% |f0| image background value (0 for black and white image with a white feature)
% |f1| image feature value (1 for black and white image with a white feature)
% Cutoff value |c| in PaLS model
% Epsilon value |eps| that is used with the smooth heaviside function
% Nu value |v| that is used with the smooth euclidean norm
% Option |opt| that determines whether centers are fixed for basis functions or whether they are allowed to float (opt=2 yields fixed centers, any other value allows them to float)

% Outputs:
% 
% Vector of values of f |vect| on the grid
function vect = f_vect_grid(p,X,Y,f0,f1,c,eps,v,opt)

%First find the sum of all basis functions on the grid
phi = phi_sum_grid(p,X,Y,v,opt);

%Use the above phi and the cutoff value to find f on the grid
z = f1*smooth_heaviside_grid(phi-c,eps) + f0*(1-smooth_heaviside_grid(phi-c,eps));

%Vectorize
vect = z(:);