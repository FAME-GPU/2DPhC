# 2DPhC
 Based on the finite element method, we developed an efficient and unified numerical algorithm for band structure calculations of 2D anisotropic photonic-crystal fibers. The main merits of this algorithm are three-fold. First, we work within the 2D primitive cell throughout, which shows that the orthogonality of the mesh is unnecessary in the discretization. Second, the quasi-periodic condition basis function involved in the FEM-based Yee’s scheme is constructed in a truly uniform manner for whatever angle between a1 and a2 of a 2D Bravais lattice. Then, using such FEM-based Yee’s scheme, the discretization turns out to have the same form for any 2D lattice. Third, based on the aid of the factorization of the stiffness matrix, we derive the nullspace free standard eigenvalue problem (SEP). It is shown that it is more efficient to solve the nullspace free SEP than the SEP with the nullspace explicitly deflated.

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
