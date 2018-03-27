function H = Heaviside(u,epsilon)
% H = 0.5*(1+2/pi*atan(u./epsilon));
x=u;
H=0.5*(erf(x/(epsilon*sqrt(2)))+1);%ÐÂµÄheavisideº¯Êý