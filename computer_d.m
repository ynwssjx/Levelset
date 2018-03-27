function d = computer_d(I,k,sigma,b,c)
%*****************数据力项更新********************************

[row,col,dim] = size(sigma);
sigma = sqrt(sigma);
for i = 1:dim
      Ici2=(I-c(:,:,i)).^2;
      n1=(1/2)*Ici2.*conv2(1./(sigma(:,:,i).^2),k,'same');
      n2=(1/2)*conv2((b.^2./(sigma(:,:,i).^2)),k,'same');
      n3=(c(:,:,i)- I).*conv2((b./(sigma(:,:,i).^2)),k,'same');
      d(:,:,i) = conv2(log(sqrt(2*pi).*sigma(:,:,i)),k,'same')+...
        n1+n2+n3;

end 
