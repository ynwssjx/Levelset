function u=CV_Evolution(u,Img,epsilon,delta_t,L,S,P,lamda1,lamda2)

u=NeumannBoundCond(u);
C=curvature_central(u); 

penalizeTerm=penalize(u);%双谷距离惩罚项函数

H_u = 0.5*(1+(2/pi)*atan(u/epsilon));
delta_H = (1/pi)*epsilon./(epsilon^2+u.^2);
c1=sum(sum(H_u.*Img))/sum(sum(H_u));
c2=sum(sum((1-H_u).*Img))/sum(sum(1-H_u));

CV_dataForce=delta_H.*(lamda1*(Img-c1).^2-lamda2*(Img-c2).^2);

% penalizeTerm=P*(4*del2(u)-C);
penalizeTerm=P*penalizeTerm;%双谷时的距离惩罚项

lengthTerm=L*delta_H.*C;
areaTerm=delta_H.*S;

u=u+delta_t*(lengthTerm-areaTerm-CV_dataForce+penalizeTerm);



function g = NeumannBoundCond(f)
[nrow,ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);  


function penalizeTerm=penalize(Phi)
[phi_x,phi_y]=gradient(Phi);
    %s=sqrt(phi_x^2+phi_y^2+1e-10);%错误写法
 s=sqrt(phi_x.*phi_x+phi_y.*phi_y);
    %*****分段函数构造*******************
    
    
 a=(s>=0)&(s<=1);
 b=(s>1);
    %*****************下面为对应的4种PS*******************************************
ps=(s-1).*b+a.*(sin(2*pi.*s)./(2*pi));%（1）原始扩散函数
    
%       ps=(s-1./((s==0)+(s~=0).*s)).*b+a.*(sin(2*pi.*s)./(2*pi));%（2）大于1 时候的扩散函数 改变 ；而小于1时 跟原始一样
%     
%      ps=(s-1).*b+(2*s.^3-3*s.^2+s).*a;% （3）大于1 时候 没有改变； 而小于1时候 改变
%     
%  ps=(s-1./((s==0)+(s~=0).*s)).*b+a.*(2*s.^3-3*s.^2+s);% （4）大于1  和  小于1  时候均改变
%       ps=(2*s.^3-3*s.^2+s).*b+a.*(sin(2*pi.*s)./(2*pi));

      
   
%********************下面为扩散率函数*****************************************
 dps=((ps==0)+(ps~=0).*ps)./((s==0)+(s~=0).*s);
 
[nxx,junk]=gradient(dps.*phi_x - phi_x);  
[junk,nyy]=gradient(dps.*phi_y - phi_y);
penalizeTerm=4*del2(Phi)+nxx+nyy;

%**************************************************************************



function k = curvature_central(u)
% compute curvature for u with central difference scheme
[ux,uy] = gradient(u);
normDu = sqrt(ux.^2+uy.^2+1e-10);
Nx = ux./normDu;
Ny = uy./normDu;
[nxx,junk] = gradient(Nx);
[junk,nyy] = gradient(Ny);
k = nxx+nyy;
