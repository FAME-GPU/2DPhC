function [V, omega, iter_number, cpu_time] = eigs_2D_PhC(kx, ky, numerator, grid_nums, epsinv_array, ...
                                                        muz_array, N_wanted, opts)                                

% solve the Helmholtz eigenvalue problems for the TE or TM waves
% kx = dot(k,a1) and ky = dot(k,a2) must be provided

% wave vector k = kx * b1 + ky * b2;
% where [b1, b2] = inv([u1, u2]).'; 
% and the factor of (2*pi) has been absorbed into kx and ky
% N_wanted is the maximum eigenvalues of interest
% if SEP_flag = true, then output of V is the periodic part of the eigenvectors of 
% B^{-1/2}*A*B^{-1/2} instead of the GEP Av = \lambda Bv
% default is periodic_flag = false, and output of V is eigenvectors of the GEP Av = \lambda Bv

global iter_count
iter_count = 0;

tmp = [kx, ky];
tmp = expofix(tmp, numerator);
expikx = tmp(1); expiky = tmp(2);

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

  B_sqrt = spdiags(1.0./sqrt(muz_array(:)), 0, N, N);
  A_scaled = B_sqrt * A * B_sqrt;  clear B_sqrt;
  A_scaled = (A_scaled + A_scaled')/2.0;   % save A.mat A_scaled;   


if abs(expikx-1.0) < 1.0e-14 && abs(expiky-1.0) < 1.0e-14
    gammaflag = 1;
  trB = sum(muz_array); 
  v_null = sqrt(muz_array(:));
  B_sqrt = spdiags(1.0./v_null, 0, N, N);
  A_scaled = B_sqrt * A * B_sqrt;  clear B_sqrt;
  A_scaled = (A_scaled + A_scaled')/2.0;   % save A.mat A_scaled;   
  [R, flag, p] = chol(A_scaled,'vector');   v_null = v_null(p);
  A_truncated = A_scaled(p(1:end-1), p(1:end-1));  
  v = rand(N,1); v = v/norm(v);
  v_null = v_null/sqrt(trB);
  opts.v0 = v - (v_null' * v) * v_null;
  dA_trnt = decomposition(A_truncated,'chol');
  Afun = @(x) GinvA_fun(x, dA_trnt, v_null);
 
  opts.dim  = 2 * (N_wanted-1); %Lanczos Subspace Dimension
  %V = zeros(N-gammaflag, N_wanted-gammaflag);
 tic;
  [V, ew, flag] = eigs(Afun, N, N_wanted-1, 'lm', ...
                   'Tolerance', opts.tol,'MaxIterations', opts.maxit, 'IsFunctionSymmetric', opts.issym, 'StartVector', opts.v0, ...
                   'SubspaceDimension', opts.dim, 'IsSymmetricDefinite', true);
  cpu_time = toc;
 %fprintf('time of eigs is %e \n',cpu_time);
else
  gammaflag = 0; opts.dim  = 2 * N_wanted;
  dA = decomposition(A_scaled,'chol');
  Afun = @(x) InvA_times_vec(x, dA);
 tic;
  [V, ew, flag] = eigs(Afun, N, N_wanted, 'lm', 'Tolerance', opts.tol, 'SubspaceDimension', opts.dim,...
                  'MaxIterations', opts.maxit, 'IsFunctionSymmetric', opts.issym );
 cpu_time = toc;
end
if any(flag),  error('eigs does not converge'); end   
omega = 1.0./sqrt(diag(ew));  
iter_number = iter_count; 


if gammaflag == 1
    %V = V - (trB .\ muz_array(1:end-1)).' * V;
    %V = bsxfun(@times, V, omega.');
%     V = [V; -muz_array(1:end-1)' * V/muz_array(end)]; 
    V = [ones(N,1)/sqrt(trB), V];
    omega = [0.0; omega];
%else
%   V = bsxfun(@times, V, omega.');
end

[omega, idx] = sort(omega, 'ascend');  %from small to large eigenvalues
V = V(:, idx); 

end

function  y = InvA_times_vec(x, dA)
     global  iter_count
      y = dA\x;
     iter_count = iter_count + 1;
end

function y =  GinvA_fun(x, Asub_decomp, v)
  global  iter_count
  y = Asub_decomp\x(1:end-1,:);
  tmp = v(1:end-1)'*y;
  tmp = v*tmp;
  y = [y - tmp(1:end-1,:); tmp(end,:)];
   iter_count = iter_count + 1;
end
