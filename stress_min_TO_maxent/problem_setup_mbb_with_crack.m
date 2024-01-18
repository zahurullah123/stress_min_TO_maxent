function [coords, conn, dm, numnode, numcell, freedofs, coords_trac, W]= ...
                                    problem_setup_mbb_with_crack(dmax, rmin)

addpath('~zahurullah/Desktop/bbb/MBB_with_crack/')           
load elements.dat
conn=elements(:,2:5)';

load nodes.dat
coords=nodes(:,2:3)'; 

numnode=length(coords(1,:)); 
numcell=length(conn(1,:)); 

figure
patch('Faces',conn','Vertices',coords','LineWidth',1,...
      'FaceColor',[0,153/255,153/255],'EdgeColor','k'); 
hold on 
plot(coords(1,:),coords(2,:),'ok','markersize',4,'markerfacecolor',[1,0,0])
xlabel('x','fontsize',14); ylabel('y','fontsize',14);
axis equal 

%Find out the dirichlet dofs fully fixed 
bc_data
plot(coords(1,bc_x),coords(2,bc_x),'ok','markersize',6,'markerfacecolor',[0,0,1])


%Find out the dirichlet dofs fully fixed 
fix_dofs1=2*bc_x-1; %right edge x only

node_LB=find(coords(1,:)==0 & coords(2,:)==0);  %left botton corner 
fix_dofs2=2*node_LB-1:2*node_LB; 
fixeddofs=[fix_dofs1 fix_dofs2]; 

alldofs = 1:2*numnode;
freedofs = setdiff(alldofs,fixeddofs);



%to find out domain of influence of the node
xy=coords(:,conn(:,1)); 
r1=sqrt( (xy(1,1)-xy(1,3))^2  + (xy(2,1)-xy(2,3))^2  ); 
r2=sqrt( (xy(1,2)-xy(1,4))^2  + (xy(2,2)-xy(2,4))^2  ); 
r=max(r1, r2); 
xspac=r;
yspac=r;
dm(1,1:numnode)=dmax*(xspac*ones(1,numnode));
dm(2,1:numnode)=dmax*(yspac*ones(1,numnode));
% a=plots_nodal_support(coords,conn,dm); 


%Traction boundary (point load)
coords_trac = coords(:,nodes_trac);
plot(coords_trac(1,:),coords_trac(2,:),'ok','markersize',8,'markerfacecolor',[0,0.8,0])
xlabel('x','fontsize',14); ylabel('y','fontsize',14);



%filter for density
%coordinate of background cells centres
gs=gauss_domain(coords,numcell,conn,1); 
coords_cells(1,:) = gs(1,:); 
coords_cells(2,:) = gs(2,:); 

%to find out domain of influence of the node
dm_cells(1,1:numcell)=rmin*(xspac*ones(1,numcell));
dm_cells(2,1:numcell)=rmin*(yspac*ones(1,numcell));

W=zeros(numcell,numcell); 
for cc=1:numcell
    gpos=coords_cells(:,cc); 
    [v,L]=nodes_in_support(numcell, coords_cells, gpos, dm_cells);
    xi=coords_cells(:,v); 
    difx=abs((gpos(1,1)-xi(1,:))); 
    dify=abs((gpos(2,1)-xi(2,:))); 
    dif=sqrt(difx.^2 + dify.^2); 
    rij=dif./sqrt(dm_cells(1,v).^2 + dm_cells(2,v).^2);
    wij=(rmin-rij)./rmin; 
    W(cc,v)=wij; 
end
W=W./sum(W,2); 
W=sparse(W); 



% figure
% patch('Faces',conn','Vertices',coords','LineWidth',1,...
%       'FaceColor',[0,153/255,153/255],'EdgeColor','k'); 
% hold on 
% plot(coords(1,:),coords(2,:),'ok','markersize',4,'markerfacecolor',[1,0,0])
% xlabel('x','fontsize',14); ylabel('y','fontsize',14);
% axis equal 

