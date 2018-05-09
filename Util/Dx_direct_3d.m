function Dxtf=Dx_direct_3d(f,H)
M=H.vol_size(1); N=H.vol_size(2); R=H.vol_size(3);
f=reshape(f,M,N,R);
Dxtf=zeros(M,N,R);

Dxtf(1,:,:)=f(1,:,:);
Dxtf(2:M,:,:)=f(2:M,:,:)-f(1:(M-1),:,:);
Dxtf=Dxtf(:);