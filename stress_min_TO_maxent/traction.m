function F = traction(F, coords, numnode, dm, coords_trac)

% Traction boundary (point load)

 ty=-1/length(coords_trac(1,:));
 
for gg=coords_trac
    gpos=gg;
    [v,L]=nodes_in_support(numnode, coords, gpos, dm);
    xi=coords(:,v);
    [phi, dphix, dphiy]=Maxent_SF(xi, gpos,  dm, v);
%     F(2*v)=F(2*v) + ty*phi'; 
    F(2*v)=F(2*v) + ty*phi'; 

end 



% % Traction boundary (point load)
% for gg=coords_trac
%     gpos=gg; 
%     
%     th=atan2(gpos(2),gpos(1)); 
% %     len=0.1; 
% %     quiver(gpos(1),gpos(2),-len*cos(th),-len*sin(th),'linewidth',4)
% 
%     [v,L]=nodes_in_support(numnode, coords, gpos, dm);
%     xi=coords(:,v);
%     [phi, dphix, dphiy]=Maxent_SF(xi, gpos,  dm, v);
%     
%     f_rad=1000;
%     tx=-f_rad*cos(th); ty=-f_rad*sin(th); 
%     F(2*v-1)=F(2*v-1) + tx*phi';     
%     F(2*v)  =F(2*v)   + ty*phi';       
%     
% %     f_tan=1000;
% %     tx=-f_tan*sin(th); ty=f_tan*cos(th); 
% %     F(2*v-1)=F(2*v-1) + tx*phi';     
% %     F(2*v)  =F(2*v)   + ty*phi';       
% end 



% for gg=gs_trac
%     gpos=gg(1:2); weight=gg(3);  jac=gg(4);
%     th=atan2(gpos(2),gpos(1)); 
% %     len=0.02; 
% %     quiver(gpos(1),gpos(2),-len*sin(th),len*cos(th),'linewidth',4)
% 
%     [v,L]=nodes_in_support(numnode, coords, gpos, dm);
%     xi=coords(:,v);
%     [phi, dphix, dphiy]=Maxent_SF(xi, gpos,  dm, v);
% 
%     f_tan=1; 
%     tx=-f_tan*sin(th); ty=f_tan*cos(th);
% 
%     F(2*v-1)=F(2*v-1) + weight*jac*tx*phi'; 
%     F(2*v)  =F(2*v)   + weight*jac*ty*phi'; 
% end 
% 




%Traction boundary (point load)
% for i=1:3
%     gpos=coords_trac(:,i);
%     [v,L]=nodes_in_support(numnode, coords, gpos, dm);
%     xi=coords(:,v);
%     [phi, dphix, dphiy]=Maxent_SF(xi, gpos,  dm, v);
%     if(i==1) %bottom
%         ty=+5;
%         F(2*v)=F(2*v) + ty*phi'; 
%     elseif (i==2) %medium 
%         tx=-10*sind(45); ty=10*cosd(45);
%         F(2*v-1)=F(2*v-1) + tx*phi';     
%         F(2*v)=F(2*v) + ty*phi';         
%     else %up  
%         tx=-5;
%         F(2*v-1)=F(2*v-1) + tx*phi'; 
%     end
% end 






% %distribute load on 2 elements
% nodes_trac=find(coords(2,:)==40); 
% nodes_trac=nodes_trac(1:3);
% 
% coords_trac = coords(:,nodes_trac); 
% F=zeros(2*numnode,1);
% 
% %FEA books 
% gauss(1,1) =-.861136311594052575224;
% gauss(1,2) =-.339981043584856264803;
% gauss(1,3) = -gauss(1,2);
% gauss(1,4) = -gauss(1,1);
% gauss(2,1) =.347854845137453857373;
% gauss(2,2) =.652145154862546142627;
% gauss(2,3) = gauss(2,2);
% gauss(2,4) = gauss(2,1);
% 
% %set up gauss points along the traction boundary
% count=1; 
% for i=1:2 % 2 elements
%     xcen=(coords_trac(1,i+1)+coords_trac(1,i))/2;
%     jcob=abs((coords_trac(1,i+1)-coords_trac(1,i))/2);
%     for j=1:4
%         mark(j)=xcen-gauss(1,j)*jcob;
%         gst(1,count)=mark(j);
%         gst(2,count)=coords_trac(2,i);
%         gst(3,count)=gauss(2,j);
%         gst(4,count)=jcob;
%         count=count+1;
%     end
% end
% 
% figure
% patch('Faces',conn','Vertices',coords','LineWidth',1,...
%       'FaceColor',[0,153/255,153/255],'EdgeColor','k'); 
% hold on 
% plot(coords(1,:),coords(2,:),'ok','markersize',8,'markerfacecolor',[1,0,0])
% xlabel('x','fontsize',14); ylabel('y','fontsize',14);
% axis equal 
% plot(gst(1,:),gst(2,:),'xg','markersize',8,'linewidth',2)
% 
% %Integrate force on the traction boundary
% for gp=gst
%     gpos=gp(1:2); weight=gp(3); jcob=gp(4);
%     [v,L]=nodes_in_support(numnode, coords, gpos, dm);
%     xi=coords(:,v);
%     [phi, dphix, dphiy]=SF_and_Drivative(xi, gpos,  dm, v);
%     tx=0;
%     ty=-1;
%     clear f1 f2 
%     for k=1:L
%         f1(2*k-1)=tx*phi(k);
%         f1(2*k)=ty*phi(k);
%     end
%     f2=f1';
%     for i=1:L
%         F(2*v(i)-1:2*v(i))=F(2*v(i)-1:2*v(i))+f2(2*i-1:2*i)*weight*jcob;
%     end
% end