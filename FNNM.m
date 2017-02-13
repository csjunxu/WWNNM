function X = FNNM( Y, NSig, m, Par )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   The algorithm solves the following optimization problem
%      1/2||Y-UV'||_F_2 + lambda/2(||U||_F_2 + ||V'||_F_2)
% warning('off');
%% initialization
[U, ~, ~] = svd(full(Y),'econ');
U = U(:,1:Par.rank);
% V = V(1:Par.rank,:);
epsl = 1e-6;

seta = mean(mean(Y.^2,1));
lambda = Par.lambdac * Par.c * NSig^2 / sqrt(seta) + epsl;

%% Alternate Direction Optimization
f_curr = 0;
for i = 1:Par.maxiter
    f_prev = f_curr;
    
    % Fix U and update V
    VT = (U' * U + lambda * eye(Par.rank)) \ (U' * Y);
    % Fix V and update U
    U = (Y * VT') / (VT * VT' + lambda * eye(Par.rank));
    % energy function
    DT = norm(Y - U * VT, 'fro') ^ 2;
    %     DT = DT(:)'*DT(:);
    RT = norm(U, 'fro') ^ 2 + norm(VT, 'fro') ^ 2;
    f_curr = 0.5 * DT + lambda * RT;
    fprintf('FNNM Energy, %d th: %2.8f\n', i, f_curr);
    if (abs(f_prev - f_curr) / f_curr < 0.001)
        break;
    end
    
end

%% result
X = U * VT;
X=X + m;

end



