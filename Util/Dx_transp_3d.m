function Dxtf=Dx_transp_3d(f,H)
M=H.vol_size(1); N=H.vol_size(2); R=H.vol_size(3);
f=reshape(f,M,N,R);
Dxtf=zeros(M,N,R);

Dxtf(M,:,:)=f(M,:,:);
Dxtf(1:(M-1),:,:)=f(1:(M-1),:,:)-f(2:M,:,:);
Dxtf=Dxtf(:);