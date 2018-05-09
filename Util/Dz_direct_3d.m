function Dztf=Dz_direct_3d(f,H)
M=H.vol_size(1); N=H.vol_size(2); R=H.vol_size(3);
f=reshape(f,M,N,R);
Dztf=zeros(M,N,R);

Dztf(:,:,1)=f(:,:,1);
Dztf(:,:,2:R)=f(:,:,2:R)-f(:,:,1:(R-1));
Dztf=Dztf(:);