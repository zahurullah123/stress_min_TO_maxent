function [f0val,df0dx,fval,dfdx]=stress_minimise(x,volfrac,pl,q,p,W,F,E0,Emin,...
    numnode,numcell,freedofs,coords,conn,gs,Dmat,nQ,L_ind,g_ind,mat_dphix,mat_dphiy,dm)

x=W*x; %filtering of x with weight W

%stress sensitivity w.r.t x using meshless
E = Emin+x'.^pl*(E0-Emin);
disp('K-Assembly')
K=assemble_K(numnode,numcell,gs,Dmat,E,nQ,L_ind,g_ind,mat_dphix,mat_dphiy);
disp('Solution')
U(freedofs) = K(freedofs,freedofs)\F(freedofs);

%calculaiton of strress (S) and relaxed stress (S_hat) and relaxed von
%misses stress(MISES_hat)
[S,S_hat,MISES,MISES_hat]=stress(x,q,numcell,Dmat,E0,U,L_ind,g_ind,mat_dphix,mat_dphiy);  

dpn_dvms=(sum(MISES_hat.^p))^(1/p-1); 
pnorm=(sum(MISES_hat.^p))^(1/p);

DvmDs=zeros(numcell,3);
for i=1:numcell
    DvmDs(i,1)=1/(2*MISES_hat(i))*(2*S_hat(i,1)-S_hat(i,2));
    DvmDs(i,2)=1/(2*MISES_hat(i))*(2*S_hat(i,2)-S_hat(i,1));
    DvmDs(i,3)=3/MISES_hat(i)*S_hat(i,3);
end
beta=zeros(numcell,1);
for i=1:numcell
    beta(i)=q*(x(i))^(q-1)*MISES_hat(i)^(p-1)*DvmDs(i,:)*S(i,:)';
end
T1=dpn_dvms*beta;

gama=zeros(2*numnode,1);
gcount=0;    
for i=1:numcell
    for gg=1:nQ*nQ
        gcount=gcount+1; 

        L=L_ind(gcount); 
        v=g_ind(gcount,1:L);
        dphix=mat_dphix(gcount,1:L);
        dphiy=mat_dphiy(gcount,1:L);
        
        Bmat=zeros(3,2*L);
        Bmat(1,1:2:end)=dphix;  Bmat(2,2:2:end)=dphiy;
        Bmat(3,1:2:end)=dphiy;  Bmat(3,2:2:end)=dphix;
       
        v_ind=zeros(1,2*L); v_ind(1:2:end)=2*v-1; v_ind(2:2:end)=2*v;
        gama(v_ind)=gama(v_ind)+x(i)^q*dpn_dvms*(E0*Dmat*Bmat)'*DvmDs(i,:)'*MISES_hat(i).^(p-1);
    end
end    
lamda=zeros(2*numnode,1);
lamda(freedofs,:)=K(freedofs,freedofs)\gama(freedofs,:);


T2=zeros(numcell,1);
gcount=0; 
for i=1:numcell
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

        Kg=jac*weight*Bmat'*E0*Dmat*Bmat;
        
        v_ind=zeros(1,2*L); v_ind(1:2:end)=2*v-1; v_ind(2:2:end)=2*v;
        ug=U(v_ind);
        
        T2(i)=T2(i)-lamda(v_ind)'*pl*x(i)^(pl-1)*Kg*ug';
    end
end
pnorm_sen=T1+T2;

colormap(jet);
subplot(2,1,1); 

patch('Faces',conn','Vertices',coords','FaceVertexCData',x,...
      'FaceColor','flat','EdgeColor','none'); axis equal off; colorbar
caxis([0 1]);

subplot(2,1,2); 
% this (0.5*sign(x-0.5)+0.5)=1 when x > 0.5 otherwise 0
MISES_plot=MISES_hat.*(0.5*sign(x-0)+0.5); 
% MISES_plot=MISES_hat; 
patch('Faces',conn','Vertices',coords','FaceVertexCData',MISES_plot,...
      'FaceColor','flat','EdgeColor','none'); 
% shading interp; 
axis equal off; colorbar; caxis([0 max(MISES_plot)]); drawnow;


%calculate volume of each cell
v_cell=zeros(numcell,1); 
gcount=0; 
for i=1:numcell
    for gg=1:nQ*nQ
        gcount=gcount+1; 
        weight=gs(3,gcount);  jac=gs(4,gcount);
        v_cell(i)=v_cell(i) + jac*weight;
    end
end
dv=v_cell/sum(v_cell); 
pnorm_sen(:) = W*pnorm_sen;
dv=W*dv; 

% fval=mean(x)-volfrac; %volume constraint
fval=sum(x.*v_cell)/sum(v_cell)-volfrac;

dfdx=dv';   %sensitivity of volume constraint
f0val=pnorm;      %objective function
df0dx=pnorm_sen; %sensitivity of objective function 





% disp('hi')
% %plot disp
% Ux=U(1:2:end);  Uy=U(2:2:end); 
% displ=zeros(2*numnode,1); 
% count=1; 
% for gg=coords
%     gpos(:,1)=gg;
%     [v,L]=nodes_in_support(numnode, coords, gpos, dm);
%     xi=coords(:,v);
%     [phi, dphix, dphiy]=Maxent_SF(xi, gpos,  dm, v);
% 
%     displ(count)=phi*Ux(v)';
%     displ(count+1)=phi*Uy(v)';
%     count=count+2;
% end
% 
% figure; scale=0.005;
% x_def(1,:)=coords(1,:)+displ(1:2:end)'*scale;
% x_def(2,:)=coords(2,:)+displ(2:2:end)'*scale;
% colormap('jet')
% ux=displ(1:2:end); uy=displ(2:2:end); usum=sqrt(ux.^2 + uy.^2); 
% patch('Faces',conn(1:4,:)','Vertices',x_def','FaceVertexCData',usum,'FaceColor','interp','EdgeColor',[0 0 0])
% shading interp; axis equal; hold on; colorbar
% %     plot(x_def(1,:),x_def(2,:),'ok','markerfacecolor','k');
% xlabel('X','fontsize',14); ylabel('Y','fontsize',14);
% grid on  
