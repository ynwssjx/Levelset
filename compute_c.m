function c = compute_c(I,k,u,b,s)
%*********************更新真实信号c***********************************
[row,col,dim] = size(u);
%ks=conv2(1./s(:,:,i),k,'same');
for i = 1:dim
    KbIu = conv2(1./s(:,:,i),k,'same').*I.*u(:,:,i)-conv2(b./s(:,:,i),k,'same').*u(:,:,i);
    Kb2u = conv2(1./s(:,:,i),k,'same').*u(:,:,i);
    de   = sum(sum(Kb2u));
    de   = de+(de==0)*eps;
    c(1:row,1:col,i) = sum(sum(KbIu))/de;
end