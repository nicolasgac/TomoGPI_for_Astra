function dJ=TomoGPI_F_Criterion(H,f,g,z,ve,vxi,M,N,R,l) %H,f,g,z,ve,vxi,M,N,R,l
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         CALCUL J_reg et dJ_reg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dg=g(:)-H*f(:);

invVedg=(ve.^(-1)).*dg(:);
dJ_MC=-H'*invVedg;
dJ_MC=reshape(dJ_MC,H.vol_size);

Dz=ihaar3_GPU(z,M,N,R,l);
er_f=f(:)-Dz(:);
dJ_reg=(vxi.^(-1)).*er_f;
dJ_reg=reshape(dJ_reg,H.vol_size);

dJ=dJ_MC+dJ_reg;

clear g_estimated dg df;









