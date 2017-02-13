function  [par]=ParSet(nSig)

par.nSig      =   nSig;                                 % Variance of the noise image
par.SearchWin =   20;                                   % Non-local patch searching window                          
par.c         =   0.4*sqrt(2);                            % Constant num for the weight vector
par.Innerloop =   2;                                    % InnerLoop Num of between re-blockmatching
par.ReWeiIter =   3;
par.min       =   16;
par.max       =   32;
par.patnum    =   32;
par.tol       =   50;
C             =   0.05/nSig;
par.delta     =   C/(1+C);   

if nSig<=20
    par.patsize       =   6;                            % Patch size
    par.Iter          =   12;                            % total iter numbers
    par.lamada        =   0.65;                         % Noise estimete parameter
    par.k             =   10;
elseif nSig <= 40
    par.patsize       =   7;
    par.Iter          =   8;
    par.lamada        =   0.60; 
    par.k             =   10;
elseif nSig<=60
    par.patsize       =   8;
    par.Iter          =   9;
    par.lamada        =   0.50; 
    par.k             =   10;
elseif nSig<=75
    par.patsize       =   9;
    par.Iter          =   9;
    par.lamada        =   0.48; 
    par.k             =   8;
else
    par.patsize       =   9;
    par.Iter          =   8;
    par.lamada        =   0.48; 
    par.k             =   8;
    
end

par.step      =   floor((par.patsize)/2-1);  

% Blockmatching and perform WNNM algorithm on all the patches in the image
% is time consuming, we just perform the blockmatching and WNNM on parts of
% patches in the image (we call these patches keypatch in explanatory notes)
% par.step is the step between each keypatch, smaller step will further
% improve the denoisng result