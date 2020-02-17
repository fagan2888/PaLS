% Joshua Enxing
% Tufts University
% Done under supervision of Misha Kilmer and Eric Miller

% Initialize parameter vector for LM

% Inputs:
% 
% Transformation |M| that is applied to unknown image
% Noisy image |img|
% Grid size |n| (i.e. the image is n pixels by n pixels)
% Number of basis functions to use |num_basis_fncs|
% Two-vector of x bounds |x_bounds| for where basis functions can be placed (i.e. if image is on [-1,1]x[-1,1], use [-0.9 0.9] so not too close to edges)
% Two-vector of y bounds |y_bounds| for where basis functions can be placed
% Cutoff value |c| in PaLS model
% Nu value |v| that is used with the smooth euclidean norm

% Outputs:
% 
% Meshgrid elements |X| and |Y|
% Vector |noisy|, which is |img| vectorized
% Initial parameter vector |p_init|
function [X,Y,noisy,p_init] = image_init_params(M,img,n,num_basis_funcs,x_bounds,y_bounds,c,v)

x = linspace(-1,1,n);
y = linspace(-1,1,n);

[X,Y] = meshgrid(x,y);

noisy = img(:);

m = ceil(sqrt(num_basis_funcs));
D = ((x_bounds(2)-x_bounds(1))/m);

%Make grid of basis functions evenly spaced 
[Z,W] = meshgrid(linspace(x_bounds(1)+D/2,x_bounds(2)-D/2,m),linspace(y_bounds(1)+D/2,y_bounds(2)-D/2,m));

p_init = zeros(4*num_basis_funcs,1);

%Want initial amplitudes to be above cutoff c
a = 2*c;

%Do this so that basis functions are like a partition of unity
r = roots([4*a -15*a 20*a -10*a 0 (a-c)]);
r = r(imag(r)==0);
r = r(r < 1);
r = r(r>0);
if x_bounds(2)==x_bounds(1)
    b = sqrt(r^2-v^2)/0.1;
else
    b = 2*sqrt(r^2-v^2)/D;
end

for i=1:m
    for j=1:m
        num = (i-1)*m+j;
        if num <= num_basis_funcs
            
            %If there is only a single basis function used, initialize it
            %at the center of mass of the image
            if num_basis_funcs == 1
                tot_mass = sum(img(:));
                [ii,jj] = ndgrid(1:size(img,1),1:size(img,2));
                xc = sum(ii(:).*img(:))/(0.5*n*tot_mass) - 1;
                yc = sum(jj(:).*img(:))/(0.5*n*tot_mass) - 1;
                p_init = [a b xc yc]';
                %p_init(4*(num-1)+1:4*(num-1)+4) = [a*(-1)^(num+1)  b Z(i,j) W(i,j)]';
                
            %For more than one basis function, alternate the amplitudes positive and negative    
            else
                p_init(4*(num-1)+1:4*(num-1)+4) = [a*(-1)^(num+1) b Z(i,j) W(i,j)]';
            end
        end
    end
end


% Uncomment if you want to see image with initial parameters

%  vect = f_vect_grid(p_init,X,Y,f0,f1,c,eps,v,opt);
%  vect = vec2mat(vect,length(X));
%  vect = vect';
%  figure;
%  imshow(vect);

noisy = M*noisy;
