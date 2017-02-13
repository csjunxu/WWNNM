function [ EPat W ] = PatEstimation( NL_mat, Self_arr, Sigma_arr, CurPat, Par )

        EPat = zeros(size(CurPat));
        W    = zeros(size(CurPat));
        for  i      =   1 : length(Self_arr)                                 % For each keypatch group
            mark    =   find(NL_mat(:,i));
            Temp    =   CurPat(:, NL_mat(mark,i));                  % Non-local similar patches to the keypatch
            M_Temp  =   repmat(mean( Temp, 2 ), [1, size(mark)]);
            Temp    =   Temp-M_Temp; % Temp 为减去均值后的相似块

            E_Temp  =   LRAppro( Temp,Sigma_arr(Self_arr(i)),M_Temp,Par.lambdac*Par.c,Par.k);
            EPat(:,NL_mat(mark,i))  = EPat(:,NL_mat(mark,i))+E_Temp;      
            W(:,NL_mat(mark,i))     = W(:,NL_mat(mark,i))+ones(Par.patsize*Par.patsize,size(NL_mat(mark,i),1));
        end
end

