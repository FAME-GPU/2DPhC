function [kx_coord, ky_coord, numerator, BZvertex_string, vertex_ind_Xcoord, GMX_idx, k_point_repeatind] ...
          = kpath_square(k_grid_num)

%gxmptidx: indices of Gamma£¬X, M and Gamma point

% sampling k-points at the Brillouin zone
% Pkx = 2*pi; Pky = 2*pi;

%if mod(k_grid_num,2) == 0
%    warning('k_grid_num should be an odd number');
%    k_grid_num = k_grid_num+1;
%end
%kx1 = linspace(0, Pkx/2.0, k_grid_num);
%ky2 = linspace(0, Pky/2.0, k_grid_num);
%longestgrid = linspace((Pkx+1i*Pky)/2.0, 0, round(k_grid_num*sqrt(2.0)));
%kx3 = real(longestgrid); ky3 = imag(longestgrid);

numerator = k_grid_num * 2;
[kx1, ky2] = deal(0 : k_grid_num);
ky1 = zeros(1, k_grid_num+1);
kx2 = repmat(k_grid_num, [1, k_grid_num+1]);

[kx3, ky3] = deal(k_grid_num : -1 : 0);
kx_coord = [kx1, kx2(2:end), kx3(2:end-1)];
ky_coord = [ky1, ky2(2:end), ky3(2:end-1)];

BZvertex_string = 'GMXG';
BZvertex_string = cellstr(BZvertex_string.');
len_path_str = length(BZvertex_string);
len_path_pt  = (len_path_str-1) * k_grid_num + 1;
vertex_ind_Xcoord = 1 : k_grid_num : len_path_pt;
% indices of key points of the Brillouin zone for drawing eigenmodes
GMX_idx = [cumsum([1, k_grid_num, k_grid_num]), 1];
k_point_repeatind = [1:length(kx_coord), 1]; 
   
end