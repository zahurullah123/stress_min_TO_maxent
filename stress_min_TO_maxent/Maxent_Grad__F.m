function [Grad_F]=Maxent_Grad__F(Lam, xi, gpos, w)

dF_dLam1=0;
dF_dLam2=0; 
Z_bar=0; 
size_xi=length(xi);

for i=1:size_xi
    xi_bar=xi(1,i)-gpos(1);
    yi_bar=xi(2,i)-gpos(2);
    Z_bar=Z_bar+w(i)*exp(-Lam(1)*xi_bar-Lam(2)*yi_bar);
    dF_dLam1=dF_dLam1-w(i)*(exp(-Lam(1)*xi_bar-Lam(2)*yi_bar)*xi_bar);
    dF_dLam2=dF_dLam2-w(i)*(exp(-Lam(1)*xi_bar-Lam(2)*yi_bar)*yi_bar);
end

Grad_F=[dF_dLam1/Z_bar dF_dLam2/Z_bar]';