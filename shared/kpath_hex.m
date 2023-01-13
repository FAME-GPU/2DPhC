function [kx_coord, ky_coord, numerator, BZvertex_string, vertex_ind_Xcoord, GMK_idx, k_point_repeatind] ...
          = kpath_hex(k_grid_num)

% GMK_idx: indices of Gamma£¬K, M and Gamma point
% sampling k-points at the Brillouin zone
% numerator = 6;
%BZ_vertex.G = [0, 0]';  BZ_vertex.K = [2, 2]';  BZ_vertex.M = [3, 0]';

numerator = k_grid_num * 6;
kx1 = 0 : 3 : 3*k_grid_num;
ky1 = zeros(1,length(kx1));

kx2 = 3*k_grid_num : -1 : 2*k_grid_num;
ky2 = 0 : 2 : 2*k_grid_num;

[kx3, ky3] = deal(2*k_grid_num : -2 : 0);
kx_coord = [kx1, kx2(2:end), kx3(2:end-1)];
ky_coord = [ky1, ky2(2:end), ky3(2:end-1)];

BZvertex_string = 'GMKG';
BZvertex_string = cellstr(BZvertex_string.');
len_path_str = length(BZvertex_string);
len_path_pt  = (len_path_str-1) * k_grid_num + 1;
vertex_ind_Xcoord = 1 : k_grid_num : len_path_pt;
% indices of key points of the Brillouin zone for drawing eigenmodes
GMK_idx = [cumsum([1, k_grid_num, k_grid_num]), 1];
k_point_repeatind = [1:length(kx_coord), 1]; 
   

   
end