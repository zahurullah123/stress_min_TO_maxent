clc; clear; close all
warning off

%input parameters (rmin is used in calculation density filter)
pl=3; volfrac=0.4; q=0.8; p=10; %pl=SIMP penality, q=stress relaxation, p=p-norm aggregation parameter
E0=1; Emin=1e-9; nu=0.3; 

rmin=1.2; 
dmax=1.5; nQ=3;

disp('Pre-Processing')
%Meshless parameters

%Plane stress D matrix (this have to be multiplied with E later)
Dmat=(1/(1-nu^2))*[1 nu 0;nu 1 0; 0 0 (1-nu)/2];

[coords, conn, dm, numnode, numcell, freedofs, ...
                coords_trac, W]= problem_setup_Lbrac(dmax,rmin); 

% [coords, conn, dm, numnode, numcell, freedofs, ...
%                 coords_trac, W]= problem_setup_step_shaped(dmax,rmin); 



%To setup integration points in each cell
gs=gauss_domain(coords,numcell,conn,nQ);

F=zeros(2*numnode,1); 
F = traction(F, coords, numnode, dm, coords_trac);

%pre-calculate shape function derivatives
[g_ind, L_ind, mat_dphix, mat_dphiy]=...
                save_SF_derivatives(numcell, numnode, gs, nQ, coords, dm);

            
m=1; epsimin=0.0000001; n=numcell;
x=volfrac*ones(n,1); 
xval=x; xold1=xval; xold2=xval;

xlb=1e-3*ones(n,1); xub=1*ones(n,1);
xmin=xlb; xmax=xub;
low=xlb; upp=xub;
c=[1e4]'; d=[0]';
a0=0; a=[0]';
raa0=0.0001; raa=0.0001;
raa0eps=0.0000001; raaeps=0.0000001;
kkttol=0;
outeriter=0; maxoutit=500;
x_his=zeros(numcell,maxoutit);
if outeriter < 0.5
    [f0val,df0dx,fval,dfdx]=stress_minimise(xval,volfrac,pl,q,p,W,F,E0,Emin,...
    numnode,numcell,freedofs,coords,conn,gs,Dmat,nQ,L_ind,g_ind,...
    mat_dphix,mat_dphiy,dm);
    innerit=0;
    outvector1 = [outeriter innerit xval'];
    outvector2 = [f0val fval'];
end

kktnorm = kkttol+1;
outit = 0;
while  outit < maxoutit
    outit   = outit+1;
    outeriter = outeriter+1;
    %%%% The parameters low, upp, raa0 and raa are calculated:
    [low,upp,raa0,raa] = ...
        asymp(outeriter,n,xval,xold1,xold2,xmin,xmax,low,upp, ...
        raa0,raa,raa0eps,raaeps,df0dx,dfdx);
    
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp] = ...
        gcmmasub(m,n,outeriter,epsimin,xval,xmin,xmax,low,upp, ...
        raa0,raa,f0val,df0dx,fval,dfdx,a0,a,c,d);
    
    xold2 = xold1;
    xold1 = xval;
    xval  = xmma;
    
    [f0val,df0dx,fval,dfdx]=stress_minimise(xval,volfrac,pl,q,p,W,F,E0,Emin,...
    numnode,numcell,freedofs,coords,conn,gs,Dmat,nQ,L_ind,g_ind,...
    mat_dphix,mat_dphiy,dm);
    
    % PRINT RESULTS
%     fprintf(' It.:%5i  P-norm Stress.:%5.4f Vol.:%5.3f \n',outit,f0val, ...
%         mean(xval(:)));
    fprintf(' It.:%5i  P-norm Stress.:%5.4f Vol.:%5.3f \n',outit,f0val, ...
        fval+volfrac);
    
    %%%% The residual vector of the KKT conditions is calculated:
    [residu,kktnorm,residumax] = ...
        kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s, ...
        xmin,xmax,df0dx,fval,dfdx,a0,a,c,d);
    outvector1 = [outeriter innerit xval'];
    outvector2 = [f0val fval'];
    x_his(:,outit)=xmma;
end