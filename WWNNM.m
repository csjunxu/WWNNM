%
function  [X]  =  WWNNM( Y, C, NSig, m, Iter )
[U,SigmaY,V] =   svd(full(Y), 'econ');
PatNum        = size(Y, 2);
Temp            =   sqrt(max( diag(SigmaY) .^ 2 - PatNum * NSig(1) ^ 2, 0 ));
W_Sam = C * sqrt(PatNum) * NSig(1) ^ 2;
for i=1:Iter
    %         W_Vec    =   (C * sqrt(PatNum) * NSig .^ 2) ./ ( Temp + eps );               % Weight vector
    W_Vec  = W_Sam ./ ( Temp + eps );
    SigmaX = soft(SigmaY, diag(W_Vec));
    Temp    = diag(SigmaX);
end
X =  U * SigmaX * V' + m;
return;
