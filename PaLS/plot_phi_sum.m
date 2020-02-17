% Joshua Enxing
% Tufts University
% Done under supervision of Misha Kilmer and Eric Miller

% Plots colormap of sum of basis functions before cutoff along with
% corresponding PaLS image

% Inputs:
%
% |f0| image background value (0 for black and white image with a white feature)
% |f1| image feature value (1 for black and white image with a white feature)
% Parameter vector |p|
% Cutoff value |c| in PaLS model
% Epsilon value |eps| that is used with the smooth heaviside function
% Nu value |v| that is used with the smooth euclidean norm
% Option |opt| that determines whether centers are fixed for basis functions or whether they are allowed to float (opt=2 yields fixed centers, any other value allows them to float)

% Outputs:
%
% Plot of colormap of phi before cutoff along with PaLS image
function plot_phi_sum(f0,f1,p,c,eps,v,opt)


phii = phi_sum_grid(p,X,Y,v,opt);
fp = f_vect_grid(p,X,Y,f0,f1,c,eps,v,opt);
mat = vec2mat(fp,length(X));
mat = mat';

figure;
subplot(1,2,1);
imagesc(phii);colorbar;
title('Phi Before Cutoff (c = 0.0067)');
subplot(1,2,2);
imshow(mat);
title('PaLS Image');