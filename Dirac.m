function DrcU = Dirac(u,epsilon)
% DrcU = (epsilon/pi)./(epsilon^2+u.^2);
x=u;
DrcU=(1/(epsilon*sqrt(2*pi)))*exp(-x.^2./(2*epsilon^2));%ÐÂµÄdiracº¯Êý