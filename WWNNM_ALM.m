
function  [X] =  WWNNM_ALM( Y, C, NSig, m, Par )
% This routine solves the following weighted nuclear norm optimization problem with column weights,
% which is more general than "lrr.m"
% min |Z|_*,P + |Y-XW|_2,1
% s.t., X = Z
% inputs:
%        Y -- D*N data matrix, D is the data dimension, and N is the number
%             of image patches.
%        W -- N*N matrix of column weights

tol = 1e-8;
maxIter = 1e2;
rho = 1.1;
max_mu = 1e10;
if ~isfield(Par, 'mu')
    Par.mu = 1e-6;
end
%% Initializing optimization variables
% intialize
% W = diag(NSig(1) ./ NSig);
W = ones(1, length(NSig));
X = zeros(size(Y));
Z = zeros(size(Y));
A = zeros(size(Y));
%% Start main loop
iter = 0;
PatNum       = size(Y,2);
TempC  = C * sqrt(PatNum) * 2 * NSig(1)^2;
while iter < maxIter
    iter = iter + 1;
    % update Z, fix X and A
    temp = X + A/Par.mu;
    [U, Sigmatemp, V] =   svd(full(temp), 'econ');
    [SigmaZ, svp] = ClosedWNNM(diag(Sigmatemp), TempC, eps);
    Z =  U(:,1:svp) * diag(SigmaZ) * V(:,1:svp)';
    % update X, fix Z and A
    X = (2 * Y * diag(W) + Par.mu * Z - A) * diag(1 ./ (W.^2 + Par.mu));
    % check the convergence conditions
    stopC = max(max(abs(Z-J)));
    if display && (iter==1 || mod(iter,10)==0 || stopC<tol)
        disp(['iter ' num2str(iter) ',mu=' num2str(Par.mu,'%2.1e') ...
            ',rank=' num2str(rank(Z,1e-4*norm(Z,2))) ',stopALM=' num2str(stopC,'%2.3e')]);
    end
    if stopC < tol
        break;
    else
        % update the multiplier A, fix Z and X
        A = A + (X - Z);
        % update the parameter mu
        Par.mu = min(max_mu, Par.mu * rho);
    end
end
X =  X + m;
return;
