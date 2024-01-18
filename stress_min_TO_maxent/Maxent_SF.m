function [PI_mxt_mat  dPI_mxt_dx  dPI_mxt_dy]=Maxent_SF(xi,gpos,dm,v)

[w  dwdx dwdy]=Maxent_Weight(xi,gpos,dm,v);
[Lam]=Maxent_Lagmul(xi,gpos, w);
%% Shape Function
Z_bar=0;
size_xi=length(xi);
for k=1:size_xi
    xi_bar=xi(1,k)-gpos(1);
    yi_bar=xi(2,k)-gpos(2); 
    Z_bar=Z_bar+w(k)*exp(-Lam(1)*xi_bar-Lam(2)*yi_bar);
end

for i=1:size_xi
    xi_bar=xi(1,i)-gpos(1);
    yi_bar=xi(2,i)-gpos(2);        
    PI_mxt_mat(i)=(w(i)*exp(-Lam(1)*xi_bar-Lam(2)*yi_bar))/Z_bar;
end


%% Shape Function Derivative
dw=[dwdx; dwdy];
dw_w=zeros(2,size_xi);
for i=1:size_xi
    if w(i)~=0%w some time goves to zeros so to avoid that case 
        dw_w(:,i)=dw(:,i)/w(i);
    else
        dw_w(:,i)=[0 0]';
    end 
end


A=0; 
for i=1:size_xi
    xi_bar=xi(1,i)-gpos(1);
    yi_bar=xi(2,i)-gpos(2);
    XY_bar=[xi_bar;yi_bar]; 
    A=A+(PI_mxt_mat(i)*XY_bar)*dw_w(:,i)';
end
H=Maxent_Hessian__F(Lam, xi, gpos, w);
invH=inv(H); 

Pi_dw=zeros(2,1); 
for i=1:size_xi
    Pi_dw=Pi_dw+PI_mxt_mat(i)*dw_w(:,i);
end

for i=1:size_xi
    xi_bar=xi(1,i)-gpos(1);
    yi_bar=xi(2,i)-gpos(2);
    XY_bar=[xi_bar;yi_bar]; 
    dPI_mxt=PI_mxt_mat(i)*((XY_bar'*(invH-invH*A))'+dw_w(:,i)-Pi_dw);   
    dPI_mxt_dx(i)=dPI_mxt(1);
    dPI_mxt_dy(i)=dPI_mxt(2);
    
end

% for i=1:size_xi
%     xi_bar=xi(1,i)-gpos(1);
%     yi_bar=xi(2,i)-gpos(2);
%     XY_bar=[xi_bar;yi_bar]; 
%     DF(:,i)=(dw(:,i)/w(i))+Lam'+(XY_bar'*(invH-invH*A))';   
% end
% 
% Dphi_1=zeros(2,1);
% for i=1:size_xi
%     Dphi_1=Dphi_1+PI_mxt_mat(i)*DF(:,i);
% end
% 
% count=1;
% for i=1:size_xi
%     xi_bar=xi(1,i)-gpos(1);
%     yi_bar=xi(2,i)-gpos(2);   
% 
%     dPI_mxt=PI_mxt_mat(i)*(DF(:,count)-Dphi_1);
% 
%     dPI_mxt_dx(i)=dPI_mxt(1);
%     dPI_mxt_dy(i)=dPI_mxt(2);
%     count=count+1;
% end
