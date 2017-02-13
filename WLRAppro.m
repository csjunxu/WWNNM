function [X ] = WLRAppro( Y,sigma,c,k )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   The algorithm solves the following optimization problem
%   1/2||Y-DA||_F_2 + lambda/2(||DW1||_F_2 + ||W2A||_F_2)

%% initialization
[row, col] = size(Y);
[D,S,A] = svd(Y,'econ');
D = D(:,1:k);
A = A(1:k,:);
tempD = zeros(size(D));
tempA = zeros(size(A));
maxiter = 50;
OMG = 1.75;
epsl = 1e-5;

lambda = c*sigma^2;
lambda1 = lambda*col/k;
lambda2 = lambda*row/k;
W1 = eye(k);W2 = eye(k);

%% Alternate Direction Optimization 
% main loop
for i = 1:maxiter
    A = (D'*D+lambda2*W1'*W1)\(D'*Y);
    tempA = OMG*A+(1-OMG)*tempA;
    A = tempA;
    W2 = diag(1./sqrt(sum(A.^2,2)));
    
    D = (Y*A')/(A*A'+lambda1*W2*W2');
    tempD = OMG*D+(1-OMG)*tempD;
    D = tempD;
    W1 = diag(1./sqrt(sum(D.^2,1)));

end

%% result
X = D*A;
    
end



