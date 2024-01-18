function H=Maxent_Hessian__F(Lam, xi, gpos, w)

Z_bar=0; 
size_xi=length(xi);

for i=1:size_xi
    xi_bar=xi(1,i)-gpos(1);
    yi_bar=xi(2,i)-gpos(2);
    Z_bar=Z_bar+w(i)*(exp(-Lam(1)*xi_bar-Lam(2)*yi_bar));
end

%H(1,1)
H11_1=0;
H11_2=0; 
for i=1:size_xi
    xi_bar=xi(1,i)-gpos(1);
    yi_bar=xi(2,i)-gpos(2);
    H11_1=H11_1+w(i)*(exp(-Lam(1)*xi_bar-Lam(2)*yi_bar))*(xi_bar^2);
    H11_2=H11_2+w(i)*(exp(-Lam(1)*xi_bar-Lam(2)*yi_bar))*(xi_bar);
end
H11=(H11_1/Z_bar)-(H11_2/Z_bar)^2; 

%H(1,2)
H12_1=0;
H12_2=0; 
H12_3=0; 
for i=1:size_xi
    xi_bar=xi(1,i)-gpos(1);
    yi_bar=xi(2,i)-gpos(2);
    H12_1=H12_1+w(i)*(exp(-Lam(1)*xi_bar-Lam(2)*yi_bar))*(xi_bar*yi_bar);
    H12_2=H12_2+w(i)*(exp(-Lam(1)*xi_bar-Lam(2)*yi_bar))*(xi_bar);
    H12_3=H12_3+w(i)*(exp(-Lam(1)*xi_bar-Lam(2)*yi_bar))*(yi_bar);
end
H12=(H12_1/Z_bar)-(H12_2/Z_bar)*(H12_3/Z_bar); 

%H(2,1)
H21=H12; 

%H(2,2)
H22_1=0;
H22_2=0; 
for i=1:size_xi
    xi_bar=xi(1,i)-gpos(1);
    yi_bar=xi(2,i)-gpos(2);
    H22_1=H22_1+w(i)*(exp(-Lam(1)*xi_bar-Lam(2)*yi_bar))*(yi_bar^2);
    H22_2=H22_2+w(i)*(exp(-Lam(1)*xi_bar-Lam(2)*yi_bar))*(yi_bar);
end
H22=(H22_1/Z_bar)-(H22_2/Z_bar)^2; 

H =[H11 H12; H21 H22];

