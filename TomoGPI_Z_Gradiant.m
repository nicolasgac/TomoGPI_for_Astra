function [z_estimated_i_plus_1,critz_out]=TomoGPI_Z_Gradiant(H,g,z_estimated_i,fv,num_iter_global,niter_gradient,critz,vxi,ve,vz,M,N,R,l,z_original)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Funtion of gradient descent for coefficients z of the HHBM method
% The ASTRA toolbox is used in this function
% Author: Li Wang
% August 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z_estimated_i_n=z_estimated_i;
numiter=niter_gradient;

% figure(6);clf,imagesc(z_estimated_i_n(:,:,R/(2^(l+1))));colormap(gray);title('z reconstructed');axis('square'),axis off;drawnow
% alphaz=0.0001;
for num_iter_gradient=1:1:numiter    
    critz.num_iter=(num_iter_global-1)*niter_gradient+num_iter_gradient;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCUL DES J,J_reg et dJ, DJ_reg ....
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %niter_done=critz.num_iter;
    %[critz.J(critz.num_iter),critz.J_MC(critz.num_iter),critz.J_reg(critz.num_iter),dJ,dJ_MC,dJ_reg]=test6_zCrit_AllVect(H,z_estimated_i_n,fv,vz,vxi,M,N,R,l);
    dJ=TomoGPI_Z_Criterion(H,g,z_estimated_i_n,fv,vz,vxi,ve,M,N,R,l);
    
%     disp('[J,J_MC,J_reg]')
%     disp([critz.J(critz.num_iter),critz.J_MC(critz.num_iter),critz.J_reg(critz.num_iter)])
    
%     figure();imagesc(dJ(:,:,R/(2^(l+1))));colormap(gray);colorbar;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCUL DU PAS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    num_alphaz=sum(dJ(:).^2);
    DdJ=ihaar3_GPU(dJ,M,N,R,l);

    YxiDdJ=(vxi.^(-0.5)).*DdJ(:);
    YzdJ=(vz.^(-0.5)).*dJ(:);
    denum_alphaz=sum(YxiDdJ(:).^2)+sum(YzdJ(:).^2);

    alphaz=num_alphaz/denum_alphaz;

    %alphaz=Backtracking_z_AllVect(H,z_estimated_i_n,fv,dJ,vz,vxi,M,N,R,l,alphaz);

    critz.alpha(critz.num_iter)=alphaz;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MISE A JOUR DE z
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    direction=-dJ;
    clear dJ;
    z_estimated_i_n=reshape(z_estimated_i_n,H.vol_size);
    z_estimated_i_n=z_estimated_i_n+alphaz*direction;
    
%     err_z=abs(z_estimated_i_n(:)-z_original(:));
%     norm_z=sum(z_original(:).^2);
%     erz=sum(err_z.^2)/norm_z;
%     disp('erz')
%     disp(erz)
    
end  

z_estimated_i_plus_1=z_estimated_i_n;
critz_out=critz;


% disp('****************************')
% disp('Descente de gradient OK !!!!')
% disp('****************************')