close all;
A = zeros(2,204);
for i = 1:1
    n = 150;
    h=2/n;
    f0 = 0;
    f1 = 1;
    num_basis = 100;
    lambda = 100;
    L = 0.5;
    eta = 2;
    c = h/2;
    eps = h/4;
    max_iters = 500;
    v = h/4;
    opt = 4;
    tol = 1E-6;
    M = speye(n*n);
    pen = 5;
    
    x = linspace(-1,1,n);
    y = linspace(-1,1,n);
    
    [X,Y] = meshgrid(x,y);

    x_bounds = [-0.9 0.9];
    y_bounds = x_bounds;
    
    m = ceil(sqrt(num_basis));
    D = ((x_bounds(2)-x_bounds(1))/m);

    [Z,W] = meshgrid(linspace(x_bounds(1)+D/2,x_bounds(2)-D/2,m),linspace(y_bounds(1)+D/2,y_bounds(2)-D/2,m));
    rect = make_rect_image(n,75,75,11,21);
    
%     u = randi(m-2)+1;
%     s = randi(m-2)+1;
%     
%     x_c = round((Z(u,s)+1)*(n/2)) + randi(2)*(-1)^(randi(2));
%     y_c = round((W(u,s)+1)*(n/2)) + randi(2)*(-1)^(randi(2));
%     width = 2*(randi(13)+3)+1;
%     height = 2*(randi(13)+3)+1;
%     rect = make_rect_image(n,x_c,y_c,width,height);
%     
%     %angle = 360*rand;
%     angle = 90;
%     rect = imrotate(rect,angle);
%     e = length(rect(:,1));
%     e = floor(e/2);
%     rect = rect(e-74:e+75,e-74:e+75);
%     
%     tot_mass = sum(rect(:));
%     [ii,jj] = ndgrid(1:size(rect,1),1:size(rect,2));
%     xc = sum(ii(:).*rect(:))/tot_mass;
%     yc = sum(jj(:).*rect(:))/tot_mass;
%     
%     rect = rect';
%     %rect = imtranslate(rect,[x_c-xc,y_c-yc]);   
%     rect = imtranslate(rect,[75-xc,75-yc]); 
    
   [p_rect,p_init,X,Y,err] = levenberg_marquardt(M,rect,n,f0,f1,num_basis,[-0.9 0.9],[-0.9 0.9],lambda,c,eps,tol,max_iters,v,opt,pen);
%    a = p_rect(1:4:end);
%    b = p_rect(2:4:end);
%    p = [a;b];
%    p = p(:);
%    A(i,:) = [m*(u-1)+s angle width height p'];
   
%    alphas = vec2mat(a,6);
%    alphas = alphas';
%    
%    betas = vec2mat(b,6);
%    betas = betas';
%     figure;
%     vect = f_vect_grid(p_init,X,Y,f0,f1,c,eps,v,opt);
%     vect = vec2mat(vect,length(X));
%     rect_mat = vect';
%     subplot(1,2,1);
%     imshow(rect);
%     title('Original image');
%     subplot(1,2,2);
%     imshow(rect_mat);
%     title('Initial PaLS Guess');
%     
%     figure;
%     vect = f_vect_grid(p_rect,X,Y,f0,f1,c,eps,v,opt);
%     vect = vec2mat(vect,length(X));
%     rect_mat = vect';
%     subplot(1,2,1);
%     imshow(rect);
%     title('Original image');
%     subplot(1,2,2);
%     imshow(rect_mat);
%     title('PaLS Image');
%     
%     phi = phi_sum_grid(p_rect,X,Y,v,1);
%     figure;
%     surf(X,Y,phi);
    
%     p_mat(opt*(i-1)+1:opt*(i-1)+opt,:) = vec2mat(p_rect,opt)';
%     
%     props = regionprops(true(size(rect)), rect, 'WeightedCentroid');
%     z = props.WeightedCentroid;
%     
%     x_c = z(1)/150;
%     y_c = z(2)/150;
%     
%     for l=1:num_basis
%         d_mat(2*(i-1)+1:2*(i-1)+3,l) = [p_mat(opt*(i-1)+1,l) p_mat(opt*(i-1)+2,l) l]';
%     end
%     
%     temp_mat = d_mat(opt*(i-1)+1:opt*(i-1)+3,:)';
%     temp_mat = sortrows(temp_mat);
%     d_mat(opt*(i-1)+1:opt*(i-1)+3,:) = temp_mat';
end