function [S,Shat,MISES,MISES_hat]=stress(x,q,numcell,Dmat,E0,U,L_ind,g_ind,mat_dphix,mat_dphiy) 
                    
%calculate stress at the centre of background cell
%using 3x3 integration points for K, we can use the the one in the centre
%here its count is 5

S=zeros(numcell,3); Shat=S; 
gcount=5;
for cc=1:numcell
    L=L_ind(gcount); 
    v=g_ind(gcount,1:L);
    
    dphix=mat_dphix(gcount,1:L);
    dphiy=mat_dphiy(gcount,1:L);

    Bmat=zeros(3,2*L);
    Bmat(1,1:2:end)=dphix;  Bmat(2,2:2:end)=dphiy;
    Bmat(3,1:2:end)=dphiy;  Bmat(3,2:2:end)=dphix;
    
    v_ind=zeros(1,2*L); v_ind(1:2:end)=2*v-1; v_ind(2:2:end)=2*v; 
    U_nv=U(v_ind);
    
    S(cc,:)=E0*Dmat*Bmat*U_nv'; %S is original stress
    Shat(cc,:)=x(cc)^q*S(cc,:); %Shat is relaxed stress
    gcount=gcount+9; %this will work only for 3x3 integration points
end

MISES=sqrt(sum(S.^2,2)-S(:,1).*S(:,2)+2.*S(:,3).^2);  %sqrt(s1^2 + s2^2 -s1*s2 - 3s12^2)
MISES_hat=sqrt(sum(Shat.^2,2)-Shat(:,1).*Shat(:,2)+2.*Shat(:,3).^2);  %sqrt(s1^2 + s2^2 -s1*s2 - 3s12^2)

