function [kx_coord, ky_coord, numerator, BZvertex_string, vertex_ind_Xcoord, GMK_idx, k_point_repeatind] ...
          = kpath_monoclinic_full_real(k_grid_num)

% GMK_idx: indices of Gamma，K, M and Gamma point
% sampling k-points at the Brillouin zone
% numerator = 1;
%BZ_vertex.G = [0, 0]';  BZ_vertex.P = [x1, y1]';  BZ_vertex.M = [x2, y2]';


a = 1;
b = sqrt(5)/2;
sin_theta = (sqrt(3)+2)/2/sqrt(5);
theta = asin(sin_theta);
R = b/a;

T = [pi/a, pi/a*(1-R*cos(theta))/R/sin(theta)];
x1 = 1/2; y1 = 1/2;

N = [pi/a*( 1 - (1-R*cos(theta) )/R/sin(theta)*cot(theta) ), pi/a/R/sin(theta)];
x2 = N(1)/(2*pi/a); y2 = (N(2)-(-2*pi/a*cot(theta))*x2 )/(2*pi/b/sin(theta));

X = [0,  pi/a/R/sin(theta)];
x3 = 0; y3= 1/2;

M = [pi/a*( (1-R*cos(theta) )/R/sin(theta)*cot(theta) - 1 ), pi/a/R/sin(theta)];
x4 = M(1)/(2*pi/a); y4 = (M(2)-(-2*pi/a*cot(theta))*x4 )/(2*pi/b/sin(theta));

P = [-pi/a, pi/a*cot(theta)];
x5=-1/2; y5=0;


numerator = k_grid_num * 1;
kx1 = linspace(0, x1*k_grid_num, k_grid_num+1);
ky1 = linspace(0, y1*k_grid_num, k_grid_num+1);

kx2 = linspace(x1*k_grid_num, x2*k_grid_num, k_grid_num+1);
ky2 = linspace(y1*k_grid_num, y2*k_grid_num, k_grid_num+1);

kx3 = linspace(x2*k_grid_num, 0, k_grid_num+1);
ky3 = linspace(y2*k_grid_num, 0, k_grid_num+1);

kx4 = linspace(0, x3*k_grid_num, k_grid_num+1);
ky4 = linspace(0, y3*k_grid_num, k_grid_num+1);

kx5 = linspace(x3*k_grid_num, x4*k_grid_num, k_grid_num+1);
ky5 = linspace(y3*k_grid_num, y4*k_grid_num, k_grid_num+1);

kx6 = linspace(x4*k_grid_num, 0, k_grid_num+1);
ky6 = linspace(y4*k_grid_num, 0, k_grid_num+1);

kx7 = linspace(0, x5*k_grid_num, k_grid_num+1);
ky7 = linspace(0, y5*k_grid_num, k_grid_num+1);

kx8 = linspace(x5*k_grid_num, x4*k_grid_num, k_grid_num+1);
ky8 = linspace(y5*k_grid_num, y4*k_grid_num, k_grid_num+1);

kx_coord = [kx1, kx2(2:end), kx3(2:end), kx4(2:end), kx5(2:end), kx6(2:end), kx7(2:end),kx8(2:end)];
ky_coord = [ky1, ky2(2:end), ky3(2:end), ky4(2:end), ky5(2:end), ky6(2:end), ky7(2:end),ky8(2:end)];

BZvertex_string = 'GTNGXMGPM';
BZvertex_string = cellstr(BZvertex_string.');
len_path_str = length(BZvertex_string);
len_path_pt  = (len_path_str-1) * k_grid_num + 1;
vertex_ind_Xcoord = 1 : k_grid_num : len_path_pt;
% indices of key points of the Brillouin zone for drawing eigenmodes
GMK_idx = [1:length(kx_coord)];
k_point_repeatind = [1:length(kx_coord)]; 
end