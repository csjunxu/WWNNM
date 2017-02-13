function  [Init_Index]  =  Block_matching(X, par,Neighbor_arr,Num_arr, SelfIndex_arr,Sigma_arr)
L         =   length(Num_arr); % 有L个要处理的图像块
Init_Index   =  zeros(par.patnum,L);

m_patch=mean(X,2);
m_blocks=repmat(m_patch,[1,size(X,2)]);
Covgs_blocks=((X-m_blocks)*(X-m_blocks)');
[U S V]=svd(Covgs_blocks);
dim=10;
P=(U(:,1:dim))';
Y=P*X;

for  i  =  1 : L
    Patch = P*X(:,SelfIndex_arr(i));
    sigma = Sigma_arr(SelfIndex_arr(i));
    Neighbors = Y(:,Neighbor_arr(1:Num_arr(i),i));    
    Dist = mean((repmat(Patch,1,size(Neighbors,2))-Neighbors).^2); 
    
    min = zeros(1,size(Neighbors,2));
    max = zeros(1,size(Neighbors,2));
    min(1:par.min) = 1;%最少相似块的个数
    max(1:par.max) = 1;%最大相似块的个数
    
    [val index] = sort(Dist);
    T=par.tol+2*sigma^2;
    mark = (val<T);
    ind = mark.*max+min;
    ind1 = find(ind);
    Init_Index(1:length(ind1),i)=Neighbor_arr(index(ind1),i); %每一列中存储该参考块的相似块下标
end
