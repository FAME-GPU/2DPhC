clc; clear; close all
addpath('shared');
addpath('case_mono_triangle');

a = 1.0; %lattice constant
b = sqrt(5)/2;
sin_theta = (sqrt(3)+2)/2/sqrt(5);
theta = asin(sin_theta);
lattice_vec = a*[1, b*cos(theta); 0, b*sin(theta)];
% u1 = lattice_vec(:,1);
% u2 = lattice_vec(:,2); 
% u1, u2 are lattice translation vectors
Gramian = dot(lattice_vec(:,[1 1 2]), lattice_vec(:,[1 2 2]), 1);
Gramian = reshape(Gramian([1 2 2 3]), 2, 2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xc = 0.5; yc = 0.5;
material_handle = @(x,y) dot([x(:)-xc, y(:)-yc], [x(:)-xc, y(:)-yc]*Gramian, 2) < radius^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nx = 400;   % grid number along x
Ny = 400;   % grid number along y
N  = Nx * Ny;
grid_nums = [Nx, Ny];  clear Nx Ny;

k_grid_num = 11;  % NO of points in each Brillouin segment
N_wanted   = 10;    % the number of bands to be calculated

[eps0, mu0] = deal(1.0);
eps1 = [5,0.5i,0;-0.5i,5,0;0,0,3];
mu1  = [7,0.8i,0;-0.8i,7,0;0,0,7];

parameter.center = [xc, yc]*lattice_vec';
parameter.type = 'polygon';
parameter.a0 = 1;
parameter.r = 0.3*a;
parameter.num_side = 3; %  the number of sides of the polygon
parameter.phi = 0;
parameter.display_grid = 'off';
parameter.lattice_vec = lattice_vec;

lattice_vec_permuted = [lattice_vec(:,2), -lattice_vec(:,1)];
Gramian_permuted = [Gramian(4), -Gramian(2); -Gramian(2), Gramian(1)];
[Area, epsinv_array, muinv_array, epsz_array, muz_array, opts] = ...
prepare_2D_PhC_handle(lattice_vec_permuted, Gramian_permuted, grid_nums, ...
               eps1, eps0, mu1, mu0, parameter);
[kx_coord, ky_coord, numerator, BZvertex_string, vertex_ind_Xcoord, ...
        GXSY_idx, k_point_repeatind] = kpath_monoclinic_full_real(k_grid_num);
%GMK_idx: indices of Gamma，K, M and Gamma point

omega_te = zeros(N_wanted, length(kx_coord));
omega_tm = zeros(N_wanted, length(kx_coord));
iternum_te = zeros(length(kx_coord),1);
TE_field = zeros(N, N_wanted, length(kx_coord));
TM_field = zeros(N, N_wanted, length(kx_coord));
iternum_tm = zeros(length(kx_coord),1);
% calculate eigenvalue

for m = 1 : length(kx_coord)
  [TE_field(:,:,m), omega_te(:,m), iternum_te(m)] = eigs_2D_PhC(kx_coord(m), ky_coord(m), ...
        numerator, grid_nums, epsinv_array, muz_array, N_wanted, opts);
  [TM_field(:,:,m), omega_tm(:,m), iternum_tm(m)] = eigs_2D_PhC(kx_coord(m), ky_coord(m), ...
        numerator, grid_nums, muinv_array, epsz_array, N_wanted, opts);  
end
fprintf('eigs finishes successfully!\n');

% draw band diagram
omega_te = omega_te / (2.0*pi*Area);
omega_tm = omega_tm / (2.0*pi*Area);     

%YIG_single_cylinder
fig_band = figure(1);
set(fig_band,'name', 'Band structure(square_lattice, Single_Cylind_TE)');
BS_ax = axes(fig_band);
Plot_2D_BS(BZvertex_string, k_grid_num, k_point_repeatind, ...
           vertex_ind_Xcoord, omega_te,  BS_ax);

fig_band = figure(2);
set(fig_band, 'name', 'Band structure(square_lattice, Single_Cylind_TM)');
BS_ax = axes(fig_band);
Plot_2D_BS(BZvertex_string, k_grid_num, k_point_repeatind, ...
                           vertex_ind_Xcoord, omega_tm,  BS_ax);
% plot iter
% fig_iter = figure(3);
% set(fig_iter, 'name', 'iter TE');
% Plot_iter(BZvertex_string, k_grid_num, k_point_repeatind, ...
%                            vertex_ind_Xcoord, iternum_te',iternum_tm');
 
   Plot_iter2 (BZvertex_string, k_grid_num, k_point_repeatind, ...
                           vertex_ind_Xcoord, iternum_te', iternum_tm' )
