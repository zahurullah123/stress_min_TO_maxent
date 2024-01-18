function Lam = Maxent_Lagmul(xi, gpos, w)

K=1;
Lam=[0 0];
Lam_rec(1,:)=Lam;

N=50; %Max number of iteration

while K<=N
    Grad_F=Maxent_Grad__F(Lam, xi, gpos, w);
    H=Maxent_Hessian__F(Lam, xi, gpos, w);
    Del_Lam=-inv(H)*Grad_F; 
%     Del_Lam=-Grad_F; 
    converg_criteria=sqrt(Grad_F(1)^2+Grad_F(2)^2);
    if converg_criteria<1e-6
%         disp('Result Converged');
%         disp(K-1)
        return
    else
        Lam1=Lam+Del_Lam';
        Lam=Lam1;
        Lam_rec(K+1,:)=Lam;
        K=K+1;
    end
end

disp('Result Not Converged');
disp(Lam)


