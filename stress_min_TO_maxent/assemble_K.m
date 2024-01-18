function K=assemble_K(numnode, numcell, gs, Dmat, E, nQ, L_ind, g_ind, mat_dphix, mat_dphiy)

K=zeros(2*numnode,2*numnode);
% Dmat=(E/(1-nu^2))*[1 nu 0;nu 1 0; 0 0 (1-nu)/2];
gcount=0; 
for cc=1:numcell
    Ecell= E(cc);
    for gg=1:nQ*nQ
        gcount=gcount+1; 
        weight=gs(3,gcount);  jac=gs(4,gcount);
        
        L=L_ind(gcount); 
        v=g_ind(gcount,1:L);
        dphix=mat_dphix(gcount,1:L);
        dphiy=mat_dphiy(gcount,1:L);
        
        Bmat=zeros(3,2*L);
        Bmat(1,1:2:end)=dphix;  Bmat(2,2:2:end)=dphiy;
        Bmat(3,1:2:end)=dphiy;  Bmat(3,2:2:end)=dphix;
        
        K1=jac*weight*Bmat'*Ecell*Dmat*Bmat;
        v_ind=zeros(1,2*L); v_ind(1:2:end)=2*v-1; v_ind(2:2:end)=2*v; 
        K(v_ind,v_ind)=K(v_ind,v_ind)+K1; 
    end
end


