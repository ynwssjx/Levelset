function sigma = compute_s(I,b,k,c,u)
% *****************更新方差（simag 是方差）**************************
 [row,col,dim] = size(c);
 I2 = I.^2;
 for i = 1:dim
     bci=b+c(:,:,i);
     ku = conv2(u(:,:,i),k,'same');
     
     KuI2 = conv2(u(:,:,i).*I2,k,'same');
     Kubci2=bci.^2.*ku;
     KIubci=-2*bci.*conv2((u(:,:,i).*I),k,'same');

     %bcKuI = -2*bc.*conv2(u(:,:,i).*I,k,'same');
    % bc2Ku = bc.^2.*conv2(u(:,:,i),k,'same'); 
    % ku = conv2(u(:,:,i),k,'same');
     nu = sum(sum(KuI2+Kubci2+KIubci));
     d =  sum(sum(ku));
     d  = d + (d==0)*eps;
     sigma(1:row,1:col,i) = nu/d;
 end
    


 
 
 


 