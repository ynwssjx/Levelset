
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
c0=2;%ˮƽ��������ʼֵ
cw=roipoly;
phi=c0*2*(0.5-cw);
contour(phi,[0,0],'r');


switch(model)
    case 1                    %CVģ��
        epsilon=1;            % Heaviside������������
        L=0.1*255^2;                 
        delta_t=0.1;
        iters=210;
        S=0;
        P = .4;
        lamda1=1;
        lamda2=1; 
        iters=50;
    case 2           %LBFģ��
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
    case 3           %Liƫ��У��KVLSģ��
        iters=20;
        
%************����һ���ֲ��ˣ����а뾶ΪR��Բ������*********
       sigma2=3;%����һ������Ҫ�Ĳ��������СҪ����ʵ��ͼ���е�inhomogenity�ı�
       K=fspecial('gaussian',4*sigma2+1,sigma2);
       one=ones(row,col);
       K_one=conv2(one,K,'smae');%���ȱ��õ�λ���
%***************************************************************

%************���������������Ȩ��ϵ��***************************
       u=.1;   %�ͷ���ϵ��
       delta_t=1;%ʱ�䲽��
       v=0.001*255^2;%������ϵ��
       epsilon=1;
%**************************************************************

%************ƫ���b�͸�������intrinsicֵci******************
b=ones(row,col);
dim=2;
for i=1:dim
    c(1:row,1:col,i)=i*(max(I(:))-min(I(:)))/dim+min(I(:));
end
%*************************************************************



    case 4           %Zhangƫ��У��ģ��(VMLSģ��)
        iters=10;
        sigma =4.5;
        K_R=fspecial('gaussian',3,sigma);
        K = ones(4*sigma+1);%The constant kerenl for SVMLS
        dim = 2; % the number of segmentation regions
        b=ones(row,col);
        s=ones(row,col,dim).*(1/sqrt(2*pi));%��̬�ֲ��ķ���
        for i = 1:dim
            c(1:row,1:col,i) = i*(max(I(:))-min(I(:)))./dim + min(I(:));
        end
        epsilon = 1;
        timestep = .1;

    case 5          %The paper from signal process
        iters=100;
        sigma2=3;%����һ������Ҫ�Ĳ��������СҪ����ʵ��ͼ���е�inhomogenity�ı�
%         K=fspecial('gaussian',4*sigma2+1,sigma2);
        K=ones(4*sigma2+1);%������ó�����ģ�壬��������������ʧ
        one=ones(row,col);
        K_one=conv2(one,K,'smae');%���ȱ��õ�λ���
        uu=.04;   %�ͷ���ϵ��
        delta_t=5;%ʱ�䲽��
%       v=0.0001*255^2;%������ϵ��
        v=.5;
        epsilon=1;
        dim=2;
        lamda1=1;
        lamda2=1;%������Ҫ�ʵ����ڶ���֮��Ĺ�ϵ������ֹ��Ŀ��������ִ�������Ҫ������
        sigma=ones(row,col,dim).*(1/sqrt(2*pi));%��̬�ֲ��ķ����ʼ��

    case 6          %MREVLS by Myself
        iters=10;
        sigma1 =.5;
        K_R=fspecial('gaussian',3,sigma1);%����ˮƽ�������Ĺ��򻯣����ø�˹���򻯶��������ã�
%ģ���ѡȡԭ���ǣ�ǿ���������أ��������أ��Լ��ϳɵ��Ƿ����ʸ߷���ͼ����
%��sigma����ȡС��ͼ��Խ�ǹ⻬��ͬ�ʣ���sigmaȡ�Ľϴ�
       sigma2 =2;
       K = ones(4*sigma2+1);%����ģ�Ͳ������Ƶĳ�����ģ��
%      K=fspecial('gaussian',round(4*sigma2+1),sigma2);%��˹��ģ��
       epsilon = 1;
       timestep = .1;
       dim = 2;
       lamda1=1;
       lamda2=1;
       mu=0.4;
       v=0.001*255^2;
       log_I=I;
 %*******************%��̬�ֲ��ķ���s��ƫ����b��c�ĳ�ʼ��*****************************
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
        u = compute_u(phi,epsilon);% ����heaviside����
        c = compute_c(I,K,u,b,s);% ����c 
        b = compute_b(I,K,u,c,s);% ����ƫ����b 
        s = compute_s(I,b,K,c,u);% ���·���s
        d = computer_d(I,K,s,b,c);%���������
        [phi,H1] = Myself_evolution(phi,d,epsilon,timestep,mu,v,lamda1,lamda2);%�����ݻ�����,����Phi
    end
    
    
    if (mod(n,10)==0)
        figure(2);
        imagesc(Img);colormap(gray);axis off;hold on;
        title(num2str(n));
        contour(phi,[0,0],'r');
        hold on;
    end
end














