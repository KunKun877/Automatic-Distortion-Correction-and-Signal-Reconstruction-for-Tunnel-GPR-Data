function [ D1 ] = drr_spatial_recon(D, MASK, N, K, Niter, eps, verb, mode, a, win_t, win_x)
%DRR_SPATIAL_RECON: Spatial domain reconstruction using DMSSA.
%   This is a modification of the original drr3drecon function,
%   removing the F-XY domain processing to work directly on a 2D matrix.
%
%  IN   D:      input 2D data (with missing values, e.g., zeros)
%       MASK:   sampling mask (1 for known, 0 for missing)
%       N:      number of singular values to be preserved
%       K:      damping factor
%       Niter:  number of maximum iterations
%       eps:    tolerance for convergence
%       verb:   verbosity flag
%       mode:   mode=1 (denoising & recon), mode=0 (reconstruction only)
%       a:      weight vector for POCS update
%       win_t:  (Optional) window size for the first dimension (time/vertical)
%       win_x:  (Optional) window size for the second dimension (space/horizontal)

% --- Parameter Defaults ---
if nargin < 3, N = 1; end
if nargin < 4, K = 4; end
if nargin < 5, Niter = 30; end
if nargin < 6, eps = 1e-5; end
if nargin < 7, verb = 0; end
if nargin < 8, mode = 0; end
if nargin < 9 || isempty(a)
    if mode == 1
        a = (Niter - (1:Niter)) / (Niter - 1); % linearly decreasing
    else
        a = ones(1, Niter);
    end
end
% Ensure 'a' has the correct length if Niter is the only thing changed
if length(a) ~= Niter
    if mode == 1
        a = (Niter-(1:Niter))/(Niter-1);
    else
        a = ones(1,Niter);
    end
end

% --- Initialization ---
[nt, nx] = size(D);
mask = MASK;

% --- Hankel Window Parameters ---
% This is the modified section. It allows for user-defined window sizes.
if nargin < 10 || isempty(win_t)
    warning('Window size for vertical axis (win_t) not provided. Defaulting to floor(nt/2)+1, which may cause memory issues for large data.');
    lx = floor(nt/2) + 1;
else
    lx = win_t;
end

if nargin < 11 || isempty(win_x)
    warning('Window size for horizontal axis (win_x) not provided. Defaulting to floor(nx/2)+1, which may cause memory issues for large data.');
    ly = floor(nx/2) + 1;
else
    ly = win_x;
end

if lx > nt, error('win_t cannot be larger than the number of rows.'); end
if ly > nx, error('win_x cannot be larger than the number of columns.'); end


% --- Main Iteration Loop ---
S_obs = D;  % The observed data is the input itself
Sn_1 = S_obs; % Initial guess

for iter = 1:Niter
    if verb
        fprintf('Iteration %d / %d ...\n', iter, Niter);
    end
    
    % Step 1: Hankel transform on the current 2D estimate
    M = P_H(Sn_1, lx, ly);
    
    % Step 2: Damped Rank Reduction
    M = P_RD(M, N, K);
    
    % Step 3: Inverse Hankel transform to get the reconstructed 2D matrix
    Sn_low_rank = P_A(M, nt, nx, lx, ly);
    
    % =========================================================================
    %  MODIFICATION START: Replaced the POCS update line with a standard form
    % =========================================================================
    % Original complex line:
    % Sn = a(iter) * S_obs + (1 - a(iter)) * mask .* Sn + (1 - mask) .* Sn;

    % New robust line:
    % This line correctly applies the low-rank estimate (Sn_low_rank) to
    % the missing parts, and reinstates the original observed data (S_obs)
    % in the known parts.
    Sn = S_obs .* mask + Sn_low_rank .* (1 - mask);
    % =========================================================================
    %  MODIFICATION END
    % =========================================================================

    
    % Check for convergence
    if norm(Sn - Sn_1, 'fro') / (norm(Sn_1, 'fro') + 1e-9) < eps
        if verb
            fprintf('Converged after %d iterations.\n', iter);
        end
        break;
    end
    
    Sn_1 = Sn;
end

D1 = Sn_1;

return;
end


% =========================================================================
%  REVISED HELPER FUNCTIONS for robust Hankel/Trajectory Matrix Operations
% =========================================================================

function [dout]=P_H(din,lx,ly)
% forming trajectory matrix (also known as Hankel matrix in SSA context)
% This is a standard, robust implementation.
[nx,ny]=size(din);
lxx = nx - lx + 1; % Number of patches vertically
lyy = ny - ly + 1; % Number of patches horizontally

dout = zeros(lx*ly, lxx*lyy, 'like', din);

patch_idx = 1;
for j = 1:lyy
    for i = 1:lxx
        patch = din(i:(i+lx-1), j:(j+ly-1));
        dout(:, patch_idx) = patch(:);
        patch_idx = patch_idx + 1;
    end
end
end

function [dout]=P_RD(din,N,K)
% Rank reduction on the trajectory matrix
[U,S,V]=svd(din,'econ'); % Using 'econ' is more efficient for non-square matrices
s_diag = diag(S);

% Ensure there is a singular value at index N+1 to avoid errors
if length(s_diag) > N
    s_noise = s_diag(N+1);
    for j=1:N
        % Apply damping
        s_diag(j) = s_diag(j)*(1 - s_noise^K / (s_diag(j)^K + 1e-16));
    end
end

% Reconstruct using only the top N singular values/vectors
num_sv_to_use = min(N, length(s_diag));
dout = U(:,1:num_sv_to_use) * diag(s_diag(1:num_sv_to_use)) * (V(:,1:num_sv_to_use)');

end

function [dout]=P_A(din_hankel,nx,ny,lx,ly)
% Averaging the trajectory matrix to output the result.
% This is the robust inverse operation corresponding to the new P_H.
lxx = nx - lx + 1; % Number of patches vertically
lyy = ny - ly + 1; % Number of patches horizontally

dout = zeros(nx, ny, 'like', din_hankel);
count = zeros(nx, ny); % Use a count matrix for accurate averaging

patch_idx = 1;
for j = 1:lyy
    for i = 1:lxx
        patch = reshape(din_hankel(:, patch_idx), lx, ly);
        dout(i:(i+lx-1), j:(j+ly-1)) = dout(i:(i+lx-1), j:(j+ly-1)) + patch;
        count(i:(i+lx-1), j:(j+ly-1)) = count(i:(i+lx-1), j:(j+ly-1)) + 1;
        patch_idx = patch_idx + 1;
    end
end

% Avoid division by zero
count(count == 0) = 1;
dout = dout ./ count;
end
