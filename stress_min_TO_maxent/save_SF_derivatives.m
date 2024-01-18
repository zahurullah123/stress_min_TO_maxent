function [g_ind, L_ind, mat_dphix, mat_dphiy]=save_SF_derivatives(numcell, numnode, gs, nQ, coords, dm)

L_ind=zeros(1,numcell*nQ*nQ); 
g_ind=zeros(numcell*nQ*nQ, 10);
mat_dphix=g_ind;
mat_dphiy=g_ind; 

gcount=0; 
for cc=1:numcell
    for gg=1:nQ*nQ
        gcount=gcount+1; gpos(:,1)=gs(1:2,gcount);
        [v,L]=nodes_in_support(numnode, coords, gpos, dm);
        xi=coords(:,v);
        [phi, dphix, dphiy]=Maxent_SF(xi, gpos,  dm, v);
        
        L_ind(gcount)=L; 
        g_ind(gcount,1:L)=v; 
        mat_dphix(gcount,1:L)=dphix; 
        mat_dphiy(gcount,1:L)=dphiy; 
    end
end
