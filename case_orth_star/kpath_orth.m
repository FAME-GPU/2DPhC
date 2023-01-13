function [kx_coord, ky_coord, numerator, BZvertex_string, vertex_ind_Xcoord, GXSY_idx, k_point_repeatind] ...
          = kpath_orth(k_grid_num)

% GMK_idx: indices of Gamma，K, M and Gamma point
% sampling k-points at the Brillouin zone
% numerator = 6;
%BZ_vertex.G = [0, 0]';  BZ_vertex.K = [2, 2]';  BZ_vertex.M = [3, 0]';

numerator = k_grid_num * 2;
kx1 = 0 : 1 : 1*k_grid_num;
ky1 = zeros(1,length(kx1));

ky2 = 0 : 1 : 1*k_grid_num;
kx2 = k_grid_num*ones(1,length(ky2));


kx3 = 1*k_grid_num: -1 : 0;
ky3 = k_grid_num*ones(1,length(kx3));

ky4 = 1*k_grid_num: -1 :0;
kx4 = zeros(1,length(ky4));



kx_coord = [kx1, kx2(2:end), kx3(2:end), kx4(2:end-1)];
ky_coord = [ky1, ky2(2:end), ky3(2:end), ky4(2:end-1)];

BZvertex_string = 'GXSYG';
BZvertex_string = cellstr(BZvertex_string.');
len_path_str = length(BZvertex_string);
len_path_pt  = (len_path_str-1) * k_grid_num + 1;
vertex_ind_Xcoord = 1 : k_grid_num : len_path_pt;
% indices of key points of the Brillouin zone for drawing eigenmodes
GXSY_idx = [cumsum([1, k_grid_num, k_grid_num, k_grid_num]), 1];
k_point_repeatind = [1:length(kx_coord), 1]; 
   
end