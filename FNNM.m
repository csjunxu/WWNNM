function X = FNNM( Y, NSig, m, Par )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   The algorithm solves the following optimization problem
%      1/2||Y-UV'||_F_2 + lambda/2(||U||_F_2 + ||V'||_F_2)
% warning('off');
%% initialization
[row, col] = size(Y);
U = randn(row, Par.rank);
tempU = zeros(size(U));
VT = randn(Par.rank, col);
tempVT = zeros(size(VT));
maxiter = 6;
OMG = 1.75;
epsl = 1e-5;

seta = mean(mean(Y.^2, 1));
lambda = Par.lambdac * Par.c * NSig^2 / sqrt(seta) + epsl;
lambda1 = lambda * col / Par.rank;
lambda2 = lambda * row / Par.rank;

%% Alternate Direction Optimization
f_curr = 0;
for i = 1:maxiter
    f_prev = f_curr;
    
    % Fix U and update V
    VT = (U' * U + lambda2 * eye(Par.rank)) \ (U' * Y);
    tempVT = OMG * VT + (1 - OMG) * tempVT;
    VT = tempVT;
    % Fix V and update U
    U = (Y * VT') / (VT * VT' + lambda1 * eye(Par.rank));
    tempU = OMG * U + (1 - OMG) * tempU;
    U = tempU;
    % energy function
    DT = 0.5 * norm(Y - U * VT, 'fro') ^ 2;
    %     DT = DT(:)'*DT(:);
    RT = lambda1 * norm(U, 'fro') ^ 2 + lambda2 * norm(VT, 'fro') ^ 2;
    f_curr = DT + RT;
%     fprintf('FNNM Energy, %d th: %2.8f\n', i, f_curr);
    if (abs(f_prev - f_curr) / f_curr < 0.001)
        break;
    end
    
end

%% result
X = U * VT;
X=X + m;

end



