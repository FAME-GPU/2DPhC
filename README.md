# 2DPhC
The goal of 2DPhC is to solve the energy band of 2D anisotropic photonic-crystal fibers and get the accurate information of its energy band as soon as possible. 

# Installation
Download 2DPhC from Github.

# Operating Instructions:
(1) First, open and just run add_folder_path.m.

(2)  Then, open case_XXX_XXX.m to set your parameters.  
The following is the description of parameters:  
1. Lattice settings  
lattice_vec: Lattice translation vector;
grid_nums: The grid numbers;  
parameter.center: Barycenter coordinates; 
parameter.func_type: Scatterer shape;
parameter.r: Radius, the distance from the barycenter to the vertex  and  the parameter for the star curve when the scatterer is a circle, quilateral triangle and rounded star, respectively;

2. Other settings   
N_wanted: Set the number of bands you want;  
k_grid_num: NO of points in each Brillouin segment;
eps1: Permittivity tensor;
mu1: Permeability tensor;

(2) Finally, just run case_XXX_XXX.m to compute the band structure.  

# Support:
If you have any questions, please let us know, we are willing to help you.  
School of Mathematics, Southeast University, Nanjing, China.  
Qing Liu: qliu@seu.edu.cn  
Hao-Nan Yang: hn_yang@seu.edu.cn
