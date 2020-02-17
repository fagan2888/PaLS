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
    
    [p_rect,p_init,X,Y,err] = levenberg_marquardt(M,rect,n,f0,f1,num_basis,[-0.9 0.9],[-0.9 0.9],lambda,c,eps,tol,max_iters,v,opt,pen);

end
