
clc;
clear all;
close all;

model=5;
Img=imread('new_img.bmp');
Img=double(Img(:,:,1));
I=Img;
[row,col]=size(I);

figure(1);
imagesc(Img);colormap(gray); axis off;hold on;
c0=2;%水平集函数初始值
cw=roipoly;
phi=c0*2*(0.5-cw);
contour(phi,[0,0],'r');


switch(model)
    case 1                    %CV模型
        epsilon=1;            % Heaviside函数参数设置
        L=0.1*255^2;                 
        delta_t=0.1;
        iters=210;
        S=0;
        P = .4;
        lamda1=1;
        lamda2=1; 
        iters=50;
    case 2           %LBF模型
        iters = 100;
        numIter=1;
        lambda1 = 1.0;
        lambda2 = 1.0;
        nu = 0.003*255*255;
        timestep = 0.1;
        mu = 1;
        epsilon = 1.0;
        sigma=5;    
        K=fspecial('gaussian',round(2*sigma)*2+1,sigma); % Gaussian kernel
        KI=conv2(Img,K,'same');  
        KONE=conv2(ones(size(Img)),K,'same');
    case 3           %Li偏差校正KVLS模型
        iters=20;
        
%************产生一个局部核（文中半径为R的圆形区域）*********
       sigma2=3;%这是一个很重要的参数，其大小要根据实际图像中的inhomogenity改变
       K=fspecial('gaussian',4*sigma2+1,sigma2);
       one=ones(row,col);
       K_one=conv2(one,K,'smae');%事先备好单位卷积
%***************************************************************

%************能量函数各个项的权重系数***************************
       u=.1;   %惩罚项系数
       delta_t=1;%时间步长
       v=0.001*255^2;%长度项系数
       epsilon=1;
%**************************************************************

%************偏差函数b和各分区的intrinsic值ci******************
b=ones(row,col);
dim=2;
for i=1:dim
    c(1:row,1:col,i)=i*(max(I(:))-min(I(:)))/dim+min(I(:));
end
%*************************************************************



    case 4           %Zhang偏差校正模型(VMLS模型)
        iters=10;
        sigma =4.5;
        K_R=fspecial('gaussian',3,sigma);
        K = ones(4*sigma+1);%The constant kerenl for SVMLS
        dim = 2; % the number of segmentation regions
        b=ones(row,col);
        s=ones(row,col,dim).*(1/sqrt(2*pi));%正态分布的方差
        for i = 1:dim
            c(1:row,1:col,i) = i*(max(I(:))-min(I(:)))./dim + min(I(:));
        end
        epsilon = 1;
        timestep = .1;

    case 5          %The paper from signal process
        iters=100;
        sigma2=3;%这是一个很重要的参数，其大小要根据实际图像中的inhomogenity改变
%         K=fspecial('gaussian',4*sigma2+1,sigma2);
        K=ones(4*sigma2+1);%必须采用常量核模板，否则轮廓最终消失
        one=ones(row,col);
        K_one=conv2(one,K,'smae');%事先备好单位卷积
        uu=.04;   %惩罚项系数
        delta_t=5;%时间步长
%       v=0.0001*255^2;%长度项系数
        v=.5;
        epsilon=1;
        dim=2;
        lamda1=1;
        lamda2=1;%根据需要适当调节二者之间的关系可以阻止非目标区域出现大量不需要的轮廓
        sigma=ones(row,col,dim).*(1/sqrt(2*pi));%正态分布的方差初始化

    case 6          %MREVLS by Myself
        iters=10;
        sigma1 =.5;
        K_R=fspecial('gaussian',3,sigma1);%用以水平集函数的规则化，采用高斯规则化抖动大，少用！
%模板的选取原则是：强度异质严重，噪声严重（自己合成的那幅异质高方差图），
%则sigma尽量取小；图像越是光滑和同质，则sigma取的较大
       sigma2 =2;
       K = ones(4*sigma2+1);%用以模型参数估计的常量核模板
%      K=fspecial('gaussian',round(4*sigma2+1),sigma2);%高斯核模板
       epsilon = 1;
       timestep = .1;
       dim = 2;
       lamda1=1;
       lamda2=1;
       mu=0.4;
       v=0.001*255^2;
       log_I=I;
 %*******************%正态分布的方差s，偏差域b和c的初始化*****************************
       s=ones(row,col,dim).*(1/sqrt(2*pi));
       b=ones(row,col);
       for i = 1:dim
          c(1:row,1:col,i) = i*(max(log_I(:))-min(log_I(:)))./dim + min(log_I(:));
       end
end

for n=1:iters
    if(model==1)%CV
        phi=CV_Evolution(phi,Img,epsilon,delta_t,L,S,P,lamda1,lamda2);
    end
    
    
    if(model==2)%LBF
        phi=EVOL_LBF(phi,Img,K,KI,KONE,nu,timestep,mu,lambda1,lambda2,epsilon,numIter);
    end
    
    
    if(model==3)%KVLS
        [phi,b,c]=KVLS_Evolution(I,phi,u,v,delta_t,epsilon,K,K_one,b,c,dim);
    end
    
    
    if (model==4)%VMLS        
       u = compute_u(phi,epsilon);% the membership function, i.e., u(:,:,1)=M1,...u(:,:,4)=M4 as Eq.(11)
       c = compute_c(I,K,u,b,s);% ci in Eq.(12)
       b = compute_b(I,K,u,c,s);% b in Eq.(12)
       s = compute_s(I,b,K,c,u);% the variance of corresponding region. see the sigma in Eq.(12)
       d = computer_d(I,K,s,b,c);%the d in Eq.(10)
       phi = VMLS_evolution(phi,d,epsilon,timestep);%Level set evolution equation
       phi = conv2(phi,K_R,'same');
    end
    
    if(model==5)
        phi=SignalProcess_Evolution(I,phi,uu,v,delta_t,epsilon,K,dim,K_one,sigma,lamda1,lamda2);
    end
    
    if(model==6)
        u = compute_u(phi,epsilon);% 更新heaviside函数
        c = compute_c(I,K,u,b,s);% 更新c 
        b = compute_b(I,K,u,c,s);% 更新偏差域b 
        s = compute_s(I,b,K,c,u);% 更新方差s
        d = computer_d(I,K,s,b,c);%数据项计算
        [phi,H1] = Myself_evolution(phi,d,epsilon,timestep,mu,v,lamda1,lamda2);%调用演化方程,更新Phi
    end
    
    
    if (mod(n,10)==0)
        figure(2);
        imagesc(Img);colormap(gray);axis off;hold on;
        title(num2str(n));
        contour(phi,[0,0],'r');
        hold on;
    end
end














