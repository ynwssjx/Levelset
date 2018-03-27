clc;
clear all;
close all;
r=pascal(5);
[m,n]=size(r);
v=var(r(:));
mean=sum(sum(r))/(m*n);
s=0;
for i=1:m
    for j=1:n
        s=s+(r(i,j)-mean)^2;
    end
end
myself_var=s/(m*n);
        

%*******************3x3µÄµ¥Î»Ãþ°å***********
kenel=3*ones(3);
r_pad=padarray(r,[1 1]);
[row,col]=size(r_pad);
for i=1:row-2
    for j=1:col-2
        padarray_sub=r_pad([i:i+2],[j:j+2]);
%         pad_mean(i,j)=mean(padarray_sub);
        pad_var(i,j)=var(padarray_sub(:));
    end
end

