% Joshua Enxing
% Tufts University
% Done under supervision of Misha Kilmer and Eric Miller

% Evaluates part of Jacobian of: noisy - M*f(p) for use in LM. 

% Inputs:
% 
% Transformation |M| that is applied to unknown image
% Parameter vector |p|
% Meshgrid elements |X| and |Y|
% |f0| image background value (0 for black and white image with a white feature)
% |f1| image feature value (1 for black and white image with a white feature)
% Cutoff value |c| in PaLS model
% Epsilon value |eps| that is used with the smooth heaviside function
% Nu value |v| that is used with the smooth euclidean norm
% Option |opt| that determines whether centers are fixed for basis functions or whether they are allowed to float (opt=2 yields fixed centers, any other value allows them to float)
% Matrix of values |phi|, values of full phi on the grid

% Outputs:
%
% |J|, part of the Jacobian with respect to a certain basis function
function J = jac_grid(M,p,X,Y,f0,f1,c,eps,v,opt,phi)

J = zeros(length(X)*length(Y),opt);

%Find value of phi, phi_val for a single basis function
[phi_val,R] = phi_sum_grid(p,X,Y,eps,0);

%Derivative of H 
d_H = (f1-f0)*smooth_heaviside_deriv_grid(phi-c,eps);

%df/d(alpha)
J1= d_H.*phi_val;

C = p(3)-X;
D = p(4)-Y;
V = C.*C+D.*D;
dphi = dphi_grid(R);

%df/d(beta)
J2 = (d_H.*p(1)*p(2).*V.*dphi)./R;

if opt == 4
    %df/dx, where x is x center of basis function
    J3 = p(1)*(p(2)^2).*(d_H.*C.*dphi)./R;
    
    %df/dy, where y is y center of basis function
    J4 = p(1)*(p(2)^2).*(d_H.*D.*dphi)./R;
    
    %Make into columns of Jacobian
    J(:,1) = J1(:);
    J(:,2) = J2(:);  
    J(:,3) = J3(:);
    J(:,4) = J4(:);
else
    %When centers are fixed, only need partials with respect to alpha and
    %beta
    J(:,1) = J1(:);
    J(:,2) = J2(:);   
end
