% Joshua Enxing
% Tufts University
% Done under supervision of Misha Kilmer and Eric Miller

% Levenberg-Marquardt algorithm for use in PaLS model

% Inputs:
% 
% Transformation |M| that is applied to unknown image
% Noisy image |img|
% Grid size |n| (i.e. the image is n pixels by n pixels)
% |f0| image background value (0 for black and white image with a white feature)
% |f1| image feature value (1 for black and white image with a white feature)
% Number of basis functions to use |num_basis_fncs|
% Two-vector of x bounds |x_bounds| for where basis functions can be placed (i.e. if image is on [-1,1]x[-1,1], use [-0.9 0.9] so not too close to edges)
% Two-vector of y bounds |y_bounds| for where basis functions can be placed
% Lambda value in LM scheme |lambda|
% Cutoff value |c| in PaLS model
% Epsilon value |eps| that is used with the smooth heaviside function
% Tolerance |tol| (stopping criterion is when decay in relative error is below tol)
% Max number of iterations |max_iters|
% Nu value |v| that is used with the smooth euclidean norm
% Option |opt| that determines whether centers are fixed for basis functions or whether they are allowed to float (opt=2 yields fixed centers, any other value allows them to float)
% L^1 penalty |pen| on the vector alpha of amplitudes of the basis functions

% Outputs:
% 
% Final parameter vector |p|
% Initial parameter vector |p_init|
% Meshgrid elements |X| and |Y|
% Vector of error at each iteration (i.e. value of the cost function ||noisy-M*f(p)||_2^2 + pen*||alpha||_1
function [p,p_init,X,Y,err] = levenberg_marquardt(M,img,n,f0,f1,num_basis_funcs,x_bounds,y_bounds,lambda,c,eps,tol,max_iters,v,opt,pen)

%Initialize parameters based on specifications
%X,Y meshgrid elements, noisy is vectorized version of noisy image img, 
%p is vector of parameters, organized by alpha_1, beta_1, x_1, y_1,
%alpha_2, beta_2, etc.
[X,Y,noisy,p] = image_init_params(M,img,n,num_basis_funcs,x_bounds,y_bounds,c,v);
p_init = p;

%Evaluate Jacobian and f on grid for p
[J,ind] = eval_jac(M,p,X,Y,f0,f1,c,eps,v,opt,pen,tol);
fp = f_vect_grid(p,X,Y,f0,f1,c,eps,v,opt);

%Calculate initial error
e_1 = norm(noisy - M*fp)^2 + pen*norm((p(1:4:end).*p(1:4:end)+eps).^(0.25))^2;
e_last = e_1;
err = zeros(max_iters,1);
err(1) = e_last;
counter = 1;
diff = 0;
diff_last = 100;

%Iterate until max_iters reached or stopping criterion
for i=1:max_iters
    
    %Stop if relative reduction in residual is less than tol
    if abs(diff-diff_last) > tol
        
        %Solve for delta in LM 
        vec = [(noisy - M*fp)' (sqrt(pen)*(p(1:4:end).*p(1:4:end)+eps).^(0.25))']';
        delta = -((J')*J+lambda*eye(opt*length(ind)))\(J'*vec);
        
        %Handle cases differently when basis centers are fixed or not
        if opt == 2
            temp = p;
            for l = 1:length(ind)
                temp(4*(ind(l)-1)+1:4*(ind(l)-1)+2) = temp(4*(ind(l)-1)+1:4*(ind(l)-1)+2) + delta(2*(l-1)+1:2*(l-1)+2);
            end
        else
            temp = p;
            for n = 1:length(ind)
                temp(1+4*(ind(n)-1):4*ind(n))= p(1+4*(ind(n)-1):4*ind(n)) + delta(1+4*(n-1):4*n);
            end
        end
        
        %Evaluate f for new p
        fp_temp = f_vect_grid(temp,X,Y,f0,f1,c,eps,v,opt);
        
        %Find new error
        e = norm(noisy - M*(fp_temp))^2 + pen*norm((temp(1:4:end).*temp(1:4:end)+eps).^(0.25))^2;
        
        %If error has decreased, shrink lambda in LM and accept new
        %parameter value, else increase lambda in LM
        if e < e_last
            [J,ind] = eval_jac(M,temp,X,Y,f0,f1,c,eps,v,opt,pen,tol);
            lambda = lambda/10;
            e_prev = e_last;
            e_last = e;
            diff_last = diff;
            diff = abs((e_prev-e_last)/e_prev);
            p = temp;
            fp = fp_temp;
            mat = vec2mat(fp,length(X));
            mat = mat';
            imshow(mat);
            counter = counter+1;
            err(counter) = e_last;
        else
            lambda = lambda*10;
        end
    else
        break;
    end
end

%Cut vector of errors at each step to the number of iterations if the
%stopping criterion was reached
err = err(1:counter);
        