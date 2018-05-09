function [f_estimated_i_plus_1,critf_out]=TomoGPI_F_Gradiant(H,f_estimated_i,g_real,z,num_iter_global,niter_gradient,critf,ve,vxi,M,N,R,l)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function of gradient descent for object f in the HHBM method
% The ASTRA toolbox is used in this function
% Author: Li Wang
% August 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_estimated_i_n=f_estimated_i;
numiter=niter_gradient;
for num_iter_gradient=1:1:numiter
   
    critf.num_iter=(num_iter_global-1)*niter_gradient+num_iter_gradient;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATE J,J_reg et dJ, DJ_reg ....
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     [critf.J(critf.num_iter),critf.J_MC(critf.num_iter),critf.J_reg(critf.num_iter),dJ,dJ_MC,dJ_reg]=test6_fCrit_AllVect(H,f_estimated_i_n,g_real,z,ve,vxi,M,N,R,l);
    dJ=TomoGPI_F_Criterion(H,f_estimated_i_n,g_real,z,ve,vxi,M,N,R,l);

%     disp('[J,J_MC,J_reg]')
%     disp([critf.J(critf.num_iter),critf.J_MC(critf.num_iter),critf.J_reg(critf.num_iter)])
%     disp('J J_MC J_reg')
%     crit.J(crit.num_iter)
%     crit.J_MC(crit.num_iter)
%     crit.J_reg(crit.num_iter)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCUL DU PAS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    num_alpha=sum(dJ(:).^2);   
    HdJ=H*dJ(:);
    YeHdJ=(ve.^(-0.5)).*HdJ(:);
    YxidJ=(vxi.^(-0.5)).*dJ(:);
    denum_alpha=sum(YeHdJ(:).^2)+sum(YxidJ(:).^2);
    alphaf=num_alpha/denum_alpha;

%     alphaf=Backtracking_f(H,f_estimated_i_n,g_real,z,dJ,ve,vxi,M,N,R,l,alphaf);
    
    critf.alpha(critf.num_iter)=alphaf;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MISE A JOUR DE f
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    direction=-dJ;
    clear dJ;
    f_estimated_i_n=reshape(f_estimated_i_n,H.vol_size);
    f_estimated_i_n=f_estimated_i_n+alphaf*real(direction);
end  

f_estimated_i_plus_1=f_estimated_i_n;
critf_out=critf;


% disp('****************************')
% disp('Descente de gradient OK !!!!')
% disp('****************************')