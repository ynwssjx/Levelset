function u = compute_u(phi1,epsilon,phi2)
%*************heavisideº¯Êý¸üÐÂ******************************
phi1=NeumannBoundCond(phi1);
H1 = Heaviside(phi1,epsilon);
if nargin == 3;    
    phi2=NeumannBoundCond(phi2);
    H2 = Heaviside(phi2,epsilon);
    u(:,:,1) = H1.*H2;
    u(:,:,2) = H1.*(1-H2);
    u(:,:,3) = (1-H1).*H2;
    u(:,:,4) = (1-H1).*(1-H2);
else
    u(:,:,1) = H1;
    u(:,:,2) = 1-H1;
end
