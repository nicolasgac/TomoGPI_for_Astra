function Dztf=Dz_transp_3d(f,H)
M=H.vol_size(1); N=H.vol_size(2); R=H.vol_size(3);
f=reshape(f,M,N,R);
Dztf=zeros(M,N,R);

Dztf(:,:,R)=f(:,:,R);
Dztf(:,:,1:(R-1))=f(:,:,1:(R-1))-f(:,:,2:R);
Dztf=Dztf(:);