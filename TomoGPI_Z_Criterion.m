function dJ=TomoGPI_Z_Criterion(H,g,z,f,vz,vxi,ve,M,N,R,l) %H,z,f,vz,vxi,M,N,R,l
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCUL J_reg et dJ_reg with standard criterion
% J(z) = (f-Dz)'*(1/Vxi)*(f-Dz) + z'*(1/Vz)*z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dJ_MC
Dz=ihaar3_GPU(z,M,N,R,l);
er_f=f(:)-Dz(:);
invVxi_erf=(vxi.^(-1)).*er_f;
dJ_MC=-haar3_GPU(invVxi_erf,M,N,R,l);

%% dJ_reg
dJ_reg=(vz(:).^(-1)).*z(:);
dJ_reg=reshape(dJ_reg,H.vol_size);

dJ=dJ_MC+dJ_reg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCUL J_reg et dJ_reg with a adapted criterion
% J(z) = (g-HDz)'*(1/Ve)*(g-HDz) + z'*(1/Vz)*z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% dJ_MC
% Dz=ihaar3_GPU(z,M,N,R,l);
% HDz=H*Dz(:);
% er_g=g(:)-HDz(:);
% vxi=vxi(1)*ones(size(ve));
% invVxi_erg=(vxi.^(-1)).*er_g;
% Ht_invVxi_erg=H'*invVxi_erg;
% dJ_MC=-haar3_GPU(Ht_invVxi_erg,M,N,R,l);
% 
% %% dJ_reg
% dJ_reg=(vz(:).^(-1)).*z(:);
% dJ_reg=reshape(dJ_reg,H.vol_size);
% 
% dJ=dJ_MC+dJ_reg;

clear g_estimated dg df;









