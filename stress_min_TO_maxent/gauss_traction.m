function [gs_trac]=gauss_traction(coords,element_outer)

% gauss(1,1) =-0.8611363115940525752239465;
% gauss(1,2) =-0.3399810435848562648026658;
% gauss(1,3) =+0.3399810435848562648026658;
% gauss(1,4) =+0.8611363115940525752239465;
% 
% gauss(2,1) = 0.3478548451374538573730639;
% gauss(2,2) = 0.6521451548625461426269361;
% gauss(2,3) = 0.6521451548625461426269361;
% gauss(2,4) = 0.3478548451374538573730639;

gauss(1,1) =-0.5773502691896257645091488;
gauss(1,2) =+0.5773502691896257645091488;

gauss(2,1) = 1;
gauss(2,2) = 1;


lele=length(element_outer(1,:)); 
gs_trac=zeros(1,2*lele); 

countg=1; 
for i=1:lele
    xe= [coords(1,element_outer(1,i))   coords(1,element_outer(2,i))];
    ye= [coords(2,element_outer(1,i))   coords(2,element_outer(2,i))];  
    
    lx= (xe(1)-xe(2)); 
    ly= (ye(1)-ye(2)); 
    jcob=sqrt(lx^2+ly^2)/2; 
    
    cenx= (xe(1)+xe(2))/2; 
    ceny= (ye(1)+ye(2))/2; 
    
    for ii = 1:2
        gs_trac(1,countg)= cenx-gauss(1,ii)*lx/2;
        gs_trac(2,countg)= ceny-gauss(1,ii)*ly/2;
        gs_trac(3,countg)= gauss(2,ii);
        gs_trac(4,countg)= jcob;
        countg=countg+1;
    end
%     figure(1); hold on
%     plot(xe,ye,'-ok','linewidth',3)
%     plot(gs_trac(1,countg-2:countg-1),gs_trac(2,countg-2:countg-1),'xb','markersize',16)
end 



