function [Area, epsinv_array, muinv_array, epsz_array, muz_array, opts] = ...
         prepare_2D_PhC_handle(Lattice_vec_permuted, Gramian_permuted, grid_nums, eps1, eps0, mu1, ...
                        mu0, parameter)

% Lattice_vec_permuted = [a2(:), -a1(:)];
% where a1, a2 are lattice translation vectors, 
% Gramian_permuted = Lattice_vec_permuted.' * Lattice_vec_permuted
% grid_nums(1), grid_nums(2) are number of grid points along a1 and a2 directions

% mu1, mu0, eps1, eps0 are scalar or 3-by-3 matrix in the Cartesian coordinate system
% mu0, eps0 : permeability and permittivity of background, e.g. mu0=1.0 ; eps0 =1.0;
% mu1, eps1 : permeability and permittivity of scatterer,  e.g. eps1 = 15*eye(3) and
% mu1  = [14, 12.4*1i, 0; -12.4*1i, 14, 0; 0, 0, 15];
%
% material_handle is a function handle to tell whether a point is in the scatterer or background
% opts is the option of eigs

Area = abs(det(Lattice_vec_permuted));

muinv_scatterer = eps_inv_cov(mu1, Lattice_vec_permuted, Gramian_permuted);
epsinv_scatterer = eps_inv_cov(eps1, Lattice_vec_permuted, Gramian_permuted);

muinv_background = eps_inv_cov(mu0, Lattice_vec_permuted, Gramian_permuted);
epsinv_background = eps_inv_cov(eps0, Lattice_vec_permuted, Gramian_permuted);

N = grid_nums(1)*grid_nums(2); % total grid numbers in unit cell;
muinv_array  = repmat(muinv_background(:), [1, N]);
epsinv_array = repmat(epsinv_background(:), [1, N]);

muz_array  = repmat(mu0(end),  [N,1]);
epsz_array = repmat(eps0(end), [N,1]);

[Y2, X2] = meshgrid(0:grid_nums(2)-1, (0:grid_nums(1)-1).'); %vertices
X2 = X2(:); Y2 = Y2(:);

% idx2 = material_handle(X2/grid_nums(1), Y2/grid_nums(2)); clear X2 Y2;

trans_mtx  =  parameter.lattice_vec';

idx2 = Material_Locate_Handle([X2/grid_nums(1), Y2/grid_nums(2)] * trans_mtx, parameter); 

% plot
switch parameter.display_grid 
    case {'on'}
        figure(2)
        trans_xy = [X2/grid_nums(1), Y2/grid_nums(2)]* trans_mtx;

        scatter( trans_xy(idx2, 1), trans_xy(idx2, 2), 5, 'fill');
        axis equal
        
end

clear X2 Y2;

%nnz(idx2)/length(idx2)
epsinv_array(:,idx2) = repmat(epsinv_scatterer(:), [1, nnz(idx2)]);
muinv_array(:,idx2)  = repmat(muinv_scatterer(:), [1, nnz(idx2)]);
epsz_array(idx2) = eps1(end);
 muz_array(idx2) = mu1(end);

epsinv_array = epsinv_array.';
if norm(epsinv_background - eye(2)*epsinv_background(1)) < 1.e-14 && ...
   norm(epsinv_scatterer - eye(2)*epsinv_scatterer(1)) < 1.e-14
   fprintf('scalar epsilon in orthogonal lattice, and our eigensolver may be inaccurate!\n');
else
% averaging between adjacent vertices along a_1
tmp = reshape(epsinv_array(:,1), grid_nums(1), grid_nums(2));
tmp = circshift(tmp, -1, 1);
epsinv_array(:,1) = (epsinv_array(:,1) + tmp(:))/2.0;

% averaging between adjacent vertices along a_2
tmp = reshape(epsinv_array(:,4), grid_nums(1), grid_nums(2));
tmp = circshift(tmp, -1, 2);
epsinv_array(:,4) = (epsinv_array(:,4) + tmp(:))/2.0;
end


muinv_array = muinv_array.';
if norm(muinv_background - eye(2)*muinv_background(1)) < 1.e-14 && ...
   norm(muinv_scatterer - eye(2)*muinv_scatterer(1)) < 1.e-14
   fprintf('scalar mu in orthogonal lattice, and our eigensolver may be inaccurate!\n');
else
% averaging between adjacent vertices along a_1
tmp = reshape(muinv_array(:,1), grid_nums(1), grid_nums(2));
tmp = circshift(tmp, -1, 1);
muinv_array(:,1) = (muinv_array(:,1) + tmp(:))/2.0;

% averaging between adjacent vertices along a_2
tmp = reshape(muinv_array(:,4), grid_nums(1), grid_nums(2));
tmp = circshift(tmp, -1, 2);
muinv_array(:,4) = (muinv_array(:,4) + tmp(:))/2.0;
clear tmp;
end

%opts is the option of eigs
   opts.dim  = 30; %Lanczos Subspace Dimension
   opts.maxit = 200; %Lanczos MaxIterations
   opts.tol = 1.0e-9; %Lanczos Tolerance
   opts.issym  = true;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function epsinv_permuted = eps_inv_cov(eps, Lattice_vec_permuted, Gramian_permuted)
   
if isscalar(eps)
    epsinv_permuted = Gramian_permuted/eps;
elseif isequal(size(eps), [3,3]) && ishermitian(eps)
    [~,flag] = chol(eps);
    if flag ~=0, error('Indefinite permeability or permittivity!'); end
    epsinv_permuted = eps(1:2,1:2)\Lattice_vec_permuted;
    epsinv_permuted = Lattice_vec_permuted' * epsinv_permuted;
    epsinv_permuted = (epsinv_permuted + epsinv_permuted')/2.0;
    [~,flag] = chol(epsinv_permuted);
    if flag ~=0, error('Indefinite permuted permeability or permittivity!'); end
else
    error('size of eps is wrong');
end

end