function [Phi,b,c]=KVLS_Evolution(I,Phi,u,v,delta_t,epsilon,K,K_one,b,c,dim);

%**************水平集演化函数*******************************
%**************by:山金孝 2011-09-29**********************


%**************下面为与Phi相关的各个函数*****************
Phi=neuman(Phi);
Curvature=curvature(Phi);
Haviside=haviside(epsilon,Phi);
Dirac=dirac(epsilon,Phi);
%*************************************************************


%***********文中公式的计算***************************************
b=compute_b(I,b,c,K,K_one,Phi,epsilon,dim);%计算文中式（23）
[c,repeat1,repeat2]=compute_c(I,b,c,K,Phi,epsilon,dim);%计算文中式（24）,为了不重复计算，此处返回了卷积值repeat1,repeat2
e=compute_e(I,K,b,c,K_one,dim,repeat1,repeat2);%计算文中式（17）
%*************************************************************


%****************能量方程中的各项*******************************
dataterm=-Dirac.*(e(:,:,1)-e(:,:,2));
lengthterm=v*Dirac.*Curvature;
dis_regu_term=u*Regulation_term(Phi);
%**************************************************************



Phi=Phi+delta_t*(dataterm+lengthterm+dis_regu_term);



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
f=(1/2)*(1+(2/pi)*atan(x./epsilon));%传统heaviside函数
% f=0.5*(erf(x./(epsilon*sqrt(2)))+1);%新的heaviside函数
%************************************************************

%**************Dirac函数计算式******************************
function f=dirac(epsilon,Phi)
x=Phi;
f=(1/pi)*(epsilon./(x.^2+epsilon^2));%传统draic函数
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

%***********************updata b*********************************
function f=compute_b(I,b,c,K,K_one,Phi,epsilon, dim)
u(:,:,1)=haviside(epsilon,Phi);
u(:,:,2)=1-haviside(epsilon,Phi);
J1=zeros(size(I));
J2=zeros(size(I));
for i=1:dim
    J1=J1+c(:,:,i).*u(:,:,i);
    J2=J2+c(:,:,i).^2.*u(:,:,i);
end
numerator=conv2((I.*J1),K,'same');
denominator=conv2(J2,K,'same');
% numerator=K_one.*(I.*J1);
% denominator=K_one.*J2;
f=numerator./denominator;
%*********************************************************************

%*********************updata c ****************************************
function [c,repeat1,repeat2]=compute_c(I,b,c,K,Phi,epsilon,dim)
repeat1=conv2(b,K,'same');
repeat2=conv2(b.^2,K,'same');
u(:,:,1)=haviside(epsilon,Phi);
u(:,:,2)=1-haviside(epsilon,Phi);
for i=1:dim
    numerator(:,:,i)=repeat1.*I.*u(:,:,i);
    denominator(:,:,i)=repeat2.*u(:,:,i);
    c(:,:,i)=sum(sum(numerator(:,:,i)))/sum(sum(denominator(:,:,i)));
end
%***********************************************************************

%*********************dataterm计算*************************************
function f=compute_e(I,K,b,c,K_one,dim,repeat1,repeat2)
for i=1:dim
    e(:,:,i)=I.^2.*K_one-2*c(:,:,i).*I.*repeat1+c(:,:,i).^2.*repeat2;
end
f=e;
%*********************************************************************





