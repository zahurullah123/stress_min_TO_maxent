function [coords, conn, dm, numnode, numcell, freedofs, coords_trac, W]= ...
                                    problem_setup_cant(dmax, rmin)

Lb=100;   Hb=30; 
nelx=100; nely=30;  

% Lb=200;   Hb=60; 
% nelx=200; nely=60;  


[x2d,y2d]=meshgrid(0:Lb/nelx:Lb, Hb:-Hb/nely:0);

size2d=size(x2d);
count=1;
numnode=size2d(1)*size2d(2);
numcell=nelx*nely;

for i=1:size2d(2)
    for j=1:size2d(1)
        coords(1,count)= x2d(j,i);
        coords(2,count)= y2d(j,i);
        node_num(j,i)=count;
        count=count+1;
    end
end

%to find the nodes for each cell
count1=1;
count2=1;

for i=1:nelx
    for j=1:nely
        conn(1,count1)=node_num(count2);
        conn(2,count1)=node_num(count2+1);
        conn(3,count1)=node_num(count2+1+(nely+1));
        conn(4,count1)=node_num(count2+(nely+1));
        count1=count1+1;
        count2=count2+1;
    end
    count2=count2+1;
end


%Find out the dirichlet dofs
%left
nodes_left=find(coords(1,:)==0); 
fixeddofs=zeros(1, 2*length(nodes_left));
fixeddofs(1:2:end)=2*nodes_left-1; 
fixeddofs(2:2:end)=2*nodes_left; 

alldofs = 1:2*numnode;
freedofs = setdiff(alldofs,fixeddofs);


%to find out domain of influence of the node
xspac=Lb/nelx;
yspac=Hb/nely;

dm(1,1:numnode)=dmax*(xspac*ones(1,numnode));
dm(2,1:numnode)=dmax*(yspac*ones(1,numnode));



%Traction boundary (point load)
nodes_trac=find(coords(2,:)==0); 
end_node_index= nodes_trac(end-3:end); 
coords_trac = coords(:,end_node_index); 

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


figure
patch('Faces',conn','Vertices',coords','LineWidth',1,...
      'FaceColor',[0,153/255,153/255],'EdgeColor','k'); 
hold on 
plot(coords(1,:),coords(2,:),'ok','markersize',4,'markerfacecolor',[1,0,0])
xlabel('x','fontsize',14); ylabel('y','fontsize',14);
axis equal 

