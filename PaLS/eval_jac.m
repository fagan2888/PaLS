% Joshua Enxing
% Tufts University
% Done under supervision of Misha Kilmer and Eric Miller

% Evaluates Jacobian of: noisy - M*f(p) for use in LM

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
% L^1 penalty |pen| on the vector alpha of amplitudes of the basis
% functions
% Tolerance |tol| (stopping criterion is when decay in relative error is below tol)

% Outputs:
%
% Jacobian |J|
% Vector of indices |ind| where the amplitudes of basis functions are
% nonzero (used to reduce Jacobian so more computationally efficient)
function [J,ind] = eval_jac(M,p,X,Y,f0,f1,c,eps,v,opt,pen,tol)
%Evaluates sum of basis functions on entire grid for this p
%Output matrix phi
[phi,~] = phi_sum_grid(p,X,Y,eps,opt);
ind = [];
counter = 1;

%Handle cases differently for fixed centers of basis functions vs not fixed
%centers
if opt == 2
    
    J = zeros(length(X)*length(Y)+length(p)/4,length(p)/2);
    
    %Have to find jacobian for each basis function separately due to
    %structure
    for i=1:length(p)/4
        
        if abs(p((i-1)*4+1)) > tol
            
            %Fill columns in Jacobian for each basis function
            J(1:length(X)*length(Y),counter:counter+1) = jac_grid(M,p((i-1)*4+1:(i-1)*4+4),X,Y,f0,f1,c,eps,v,opt,phi);
            counter = counter + 2;
            ind = [ind;i];
        end
    end
    col = counter - 1;
else
    J = zeros(length(X)*length(Y)+length(p)/4,length(p));
    
    %Have to find jacobian for each basis function separately due to
    %structure
    for i=1:length(p)/4
        if abs(p((i-1)*4+1)) > tol
            
            %Fill columns in Jacobian for each basis function
            J(1:length(X)*length(Y),counter:counter+3) = jac_grid(M,p((i-1)*4+1:(i-1)*4+4),X,Y,f0,f1,c,eps,v,opt,phi);
            counter = counter + 4;
            ind = [ind;i];
        end
    end
    col = counter - 1;
end

J = J(:,1:col);

%Have to make part of Jacobian that comes from the L^1 penalty term
alphas = p(1:4:end);
alphas = sqrt(pen)*0.5*alphas.*(alphas.*alphas + eps).^(-0.75);
if opt == 2
    for k=1:length(ind)
        J(length(X)*length(Y)+ind(k),2*(k-1)+1) = alphas(ind(k));
    end
else
    for k=1:length(ind)
        J(length(X)*length(Y)+ind(k),4*(k-1)+1) = alphas(ind(k));
    end
end

%Since function we are finding Jacobian of is: noisy - M*f(p), and J is the
%Jacobian of f currently, we just multiply by -M
J(1:length(X)*length(Y),:) = -M*J(1:length(X)*length(Y),:);

