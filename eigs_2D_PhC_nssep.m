function [V, omega, iter_number, cpu_time] = eigs_2D_PhC_nssep(kx, ky, numerator, grid_nums, epsinv_array, ...
                                                        muz_array, N_wanted, shift_minus, opts)
                                  
global iter_count                                  
iter_count = 0;

if sign(shift_minus) ~= -1
    error('shift_minus should be a negative real number!');
end    

tmp = [kx, ky];
tmp = expofix(tmp, numerator);
expikx = tmp(1); expiky = tmp(2);

if ~(abs(expikx-1.0) < 1.0e-14 && abs(expiky-1.0) < 1.0e-14)
    error('Pls switch to eigs_2D_PhC.m!');
end

N = prod(grid_nums); % total grid numbers in unit cell;
c1col = [-1.0; zeros(grid_nums(1)-2,1); expikx];
c1row = [-1.0, 1.0, zeros(1,grid_nums(1)-2)];
K1minus = sptoeplitz(c1col,c1row)*grid_nums(1);   %C_1 matrix is I_{n_2} \otimes K1minus
% partial derivative w.r.t x

c1col = [0.0;  1.0; zeros(grid_nums(1)-3,1); -expikx];
c1row = [0.0, -1.0, zeros(1,grid_nums(1)-3), conj(expikx)];
K1s_minus_K1 = sptoeplitz(c1col,c1row)*grid_nums(1);

clear c1row c1col;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c2col = [-1.0; zeros(grid_nums(2)-2,1); expiky];
c2row = [-1.0, 1.0, zeros(1,grid_nums(2)-2)];
K2minus = sptoeplitz(c2col,c2row)*grid_nums(2);   %C_2 matrix is  K2minus \otimes I_{n_1}
% partial derivative w.r.t y

c2col = [0.0; -1.0; zeros(grid_nums(2)-3,1); expiky];
c2row = [0.0, 1.0, zeros(1,grid_nums(2)-3), -conj(expiky)];
K2_minus_K2s = sptoeplitz(c2col,c2row)*grid_nums(2); 

clear c2row c2col;

C1 = kron(speye(grid_nums(2)), K1minus);
C2 = kron(K2minus, speye(grid_nums(1)));
CC1 = kron(speye(grid_nums(2)), K1s_minus_K1); 
CC2s = kron(K2_minus_K2s, speye(grid_nums(1)));
clear K1minus K2minus K1s_minus_K1 K2_minus_K2s; 

A1 = C1' * spdiags(epsinv_array(:,1), 0, N, N);  A1 = A1 * C1; %diagonal part
A2 = C2' * spdiags(epsinv_array(:,4), 0, N, N);  A2 = A2 * C2; %diagonal part
A  = sparse(A1 + A2); clear A1 A2;
A  = (A + A')/2.0;

if any(abs(epsinv_array(:,3)) > 1.0e-14)
   A3 = CC1 * spdiags(epsinv_array(:,3), 0, N, N); A3 = A3 * CC2s; %off-diagonal part
   A  = A + (A3 + A3')/4.0;  clear A3;
end

%save A.mat A;   
%save B.mat muz_array;

  B_sqrt = spdiags(1.0./sqrt(muz_array(:)), 0, N, N);
  A_scaled = B_sqrt * A * B_sqrt;  clear B_sqrt;
  A_scaled = (A_scaled + A_scaled')/2.0;   % save A.mat A_scaled;   
  dA = decomposition(A_scaled-shift_minus*speye(N),'chol');
  Afun = @(x) InvA_times_vec(x, dA);
  opts.dim  = 2*N_wanted; %Lanczos Subspace Dimension

tic;
[V, ew, flag] = eigs(Afun, N, N_wanted, 'lm', ...
                'Tolerance', opts.tol, 'MaxIterations', opts.maxit, 'IsFunctionSymmetric', opts.issym, ...
                'SubspaceDimension', opts.dim, 'IsSymmetricDefinite', true);
cpu_time = toc;
%fprintf('time of eigs is %e \n',cpu_time);
iter_number = iter_count;

ew = 1./diag(ew)+shift_minus;
[ew, idx] = sort(ew, 'ascend');  %from small to large eigenvalues
V = V(:, idx); 
omega = sqrt(abs(ew));

end

function  y = InvA_times_vec(x, dA)
     global  iter_count
     y = dA\x;
     iter_count = iter_count + 1;
end
