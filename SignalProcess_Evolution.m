function Phi=SignalProcess_Evolution(I,Phi,uu,v,delta_t,epsilon,K,dim,K_one,sigma,lamda1,lamda2)

%**************水平集演化函数*******************************
%**************by:山金孝 2011-11-22**********************


%**************下面为与Phi相关的各个函数*****************
Phi=neuman(Phi);
Curvature=curvature(Phi);
Haviside=haviside(epsilon,Phi);
Dirac=dirac(epsilon,Phi);
%*************************************************************


%***********文中公式的计算***************************************
u=compute_u(I,K,Phi,dim,epsilon,K_one,sigma);%**************计算文中公式（22）
sigma=compute_sigma(u,I,K,Phi,dim,epsilon);%************计算文中公式（23）
e=compute_e(I,u,sigma,K,dim);%*****************计算（25）和（26）
%*************************************************************


%****************能量方程中的各项*******************************
dataterm=-Dirac.*(lamda1*e(:,:,1)-lamda2*e(:,:,2));
lengthterm=v*Dirac.*Curvature;
% dis_regu_term=uu*Regulation_term(Phi);
dis_regu_term=uu*(4*del2(Phi)-Curvature);
%**************************************************************



Phi=Phi+delta_t*(dataterm+lengthterm+dis_regu_term);
%  Phi=Phi+delta_t*(dataterm+lengthterm);


%**********************纽曼边界（自然边界）********************
function f=neuman(Phi)
p=Phi;
[r,c]=size(p);
p([1,r],[1,c])=Phi([3,end-2],[3,c-2]);
p([1,r],[2:end-1])=Phi([3,end-2],[2:end-1]);
p([2:end-1],[1,c])=Phi([2:end-1],[3,c-2]);
f=p;
%************************************************************

%************************曲率计算函数**************************
function f=curvature(Phi)
[nx,ny]=gradient(Phi);
s=sqrt(nx.^2+ny.^2);
s=(s==0).*eps+s;
Nx=nx./s;
Ny=ny./s;
[Nxx,Nxy]=gradient(Nx);
[Nyx,Nyy]=gradient(Ny);
f=Nxx+Nyy;
%**************************************************************

%****************Hevaside函数计算式****************************
function f=haviside(epsilon,Phi)
x=Phi;
f=(1/2)*(1+(2/pi)*atan(x./epsilon));
% f=0.5*(erf(x./(epsilon*sqrt(2)))+1);%新的heaviside函数
%************************************************************

%**************Dirac函数计算式******************************
function f=dirac(epsilon,Phi)
x=Phi;
f=(1/pi)*(epsilon./(x.^2+epsilon^2));
% f=(1/(epsilon*sqrt(2*pi)))*exp(-x.^2./(2*epsilon^2));%新的dirac函数
%*************************************************************

%***************距离规则项的计算******************************
function f=Regulation_term(Phi)
[nx,ny]=gradient(Phi);
s=sqrt(nx.^2+ny.^2);
a=(s>=0)&(s<=1);
b=(s>1);
ps=(s-1./((s==0)+(s~=0).*s)).*b+a.*(2*s.^3-3*s.^2+s);
dps=((ps==0)+(ps~=0).*ps)./((s==0)+(s~=0).*s);
[nxx,nxy]=gradient(dps.*nx-nx);
[nyx,nyy]=gradient(dps.*ny-ny);
f=4*del2(Phi)+nxx+nyy;
%***************************************************************

%***********************updata u*****************************************
%************************注意文中计算u的公式是错误的********************
function f=compute_u(I,k,Phi,dim,epsilon,K_one,sigma)
[row,col]=size(I);
m(1:row,1:col,1)=haviside(epsilon,Phi);
m(1:row,1:col,2)=1-haviside(epsilon,Phi);

for i=1:dim
    sigma(1:row,1:col,i)=sigma(1:row,1:col,i)+(sigma(1:row,1:col,i)==0)*eps;
    numerator(1:row,1:col,i)=conv2((I.*m(1:row,1:col,i))./sigma(1:row,1:col,i),k,'same');
    denominator(1:row,1:col,i)=conv2(m(1:row,1:col,i)./sigma(1:row,1:col,i),k,'same');
    denominator(1:row,1:col,i)= denominator(1:row,1:col,i)+(denominator(1:row,1:col,i)==0)*eps;
    %*************u是变量**************************
    u(1:row,1:col,i)=numerator(1:row,1:col,i)./denominator(1:row,1:col,i);
   
    %*************u是常量**************************************8
%    de=sum(sum(denominator(1:row,1:col,i)));
%    de=de+(de==0)*eps;
%    u(1:row,1:col,i)=sum(sum(numerator(1:row,1:col,i)))/de;
    %*******************************************
end

% for i=1:dim
%     numerator=sum(sum(I.*m(:,:,i).*K_one));
%     denominator=sum(sum(m(:,:,i).*K_one));
%     denominator=denominator+(denominator==0)*eps;
%     u(1:row,1:col,i)=numerator/denominator;
% end

f=u;

%*********************************************************************

%*********************updata sigma ****************************************
%**********************注意文中计算sigma的公式是错误的，少一个积分符号**********************
function f=compute_sigma(u,I,K,Phi,dim,epsilon)
[row,col]=size(I);
I2=I.^2;
k=K;
m(1:row,1:col,1)=haviside(epsilon,Phi);
m(1:row,1:col,2)=1-haviside(epsilon,Phi);
for i=1:dim
    first(1:row,1:col,i)=u(1:row,1:col,i).^2.*conv2(m(1:row,1:col,i),k,'same');
    second(1:row,1:col,i)=conv2((I2.*m(1:row,1:col,i)),k,'same');
    third(1:row,1:col,i)=-2*u(1:row,1:col,i).*conv2((I.*m(1:row,1:col,i)),k,'same');
    denominator(1:row,1:col,i)=conv2(m(1:row,1:col,i),k,'same');
    denominator(1:row,1:col,i)=denominator(1:row,1:col,i)+(denominator(1:row,1:col,i)==0)*eps;
    
      %******sigma是变量******************************
    numerator(1:row,1:col,i)=first(1:row,1:col,i)+second(1:row,1:col,i)+third(1:row,1:col,i);
    sigma(1:row,1:col,i)=numerator(1:row,1:col,i)./denominator(1:row,1:col,i);
    %********************************************************************
    
      %******sigma是常量******************************
%     numerator=sum(sum(first(1:row,1:col,i)+second(1:row,1:col,i)+third(1:row,1:col,i)));
%     de=sum(sum(denominator(1:row,1:col,i)));
%     de=de+(de==0)*eps;
%     sigma(1:row,1:col,i)=numerator/de;
     %*****************************************************************
end
f=sigma;

%***********************************************************************

%*********************dataterm计算*************************************
%*********************计算e的公式也是错误的，文中公式将“sqrt(2*pi)”吞掉了！
function f=compute_e(I,u,sigma,K,dim)
[row,col]=size(I);
k=K;
sigma=sqrt(sigma);
sigma=sigma+(sigma==0)*eps;
for i=1:dim
    first(1:row,1:col,i)=conv2(log(sqrt(2*pi).*sigma(1:row,1:col,i)),k,'same');
    second(1:row,1:col,i)=(1/2)*conv2((u(1:row,1:col,i).^2./sigma(1:row,1:col,i).^2),k,'same');
    third(1:row,1:col,i)=(1/2)*I.^2.*conv2(1./sigma(1:row,1:col,i).^2,k,'same');
    fourth(1:row,1:col,i)=-I.*conv2(u(1:row,1:col,i)./sigma(1:row,1:col,i).^2,k,'same');
    e(1:row,1:col,i)=first(1:row,1:col,i)+second(1:row,1:col,i)+third(1:row,1:col,i)+fourth(1:row,1:col,i);
end
f=e;

%*********************************************************************





