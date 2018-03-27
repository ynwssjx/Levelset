function b = compute_b(I,k,u,c,sigma)
%*************¸üÐÂÆ«²îÓò***********************************************

[row,col,dim] = size(u);

J1 = zeros(row,col);
J2 = zeros(row,col);

for i = 1:dim
    
    reverse_s = 1./sigma(:,:,i);
    ku=conv2(u(:,:,i),k,'same');
   
    J1 = J1 + (conv2(I.*u(:,:,i),k,'same')-c(:,:,i).*ku).*reverse_s;
    J2 = J2 + ku.*reverse_s;
end

b = J1./(J2+(J2==0).*eps);
