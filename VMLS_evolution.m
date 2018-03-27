function phi1 = VMLS_evolution(phi1,d,epsilon,timestep)

    phi1=NeumannBoundCond(phi1);
    H1 = Heaviside(phi1,epsilon);
    Dirac1 = Dirac(phi1,epsilon);
 
%     phi2=NeumannBoundCond(phi2);
%     H2 = Heaviside(phi2,epsilon);
%     Dirac2 = Dirac(phi2,epsilon);
%     DataF1 = -(d(:,:,1)-d(:,:,2)-d(:,:,3)+d(:,:,4)).*H2-(d(:,:,2)-d(:,:,4));
     DataF=d(:,:,2)-d(:,:,1);

    phi1 = phi1 +timestep*(DataF.*Dirac1); 
%     DataF2 = -(d(:,:,1)-d(:,:,2)-d(:,:,3)+d(:,:,4)).*H1-(d(:,:,3)-d(:,:,4));
%     phi2 = phi2 +timestep*(DataF2.*Dirac2);
