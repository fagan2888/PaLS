README.txt

levenbergmarquardt.m is the main m-file for optimization, takes inputs listed in file

The following m-files are used in the execution of levenbergmarquardt.m:	
	
phi_sum_grid.m, jac_grid.m, image_init_params.m, f_vect_grid.m, eval_jac.m, dphi_grid.m, smooth_heaviside_grid.m,
smooth_heaviside_deriv_grid.m, smooth_euclidean_dx1_grid.m, smooth_euclidean_dx2_grid.m,
smooth_abs.m, and smooth_abs_deriv.m all

The remaining m-files are stand-alone and are used to produce images:

plot_phi.m, plot_phi_sum.m, make_circle.m, make_rect_image.m