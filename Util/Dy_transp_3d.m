function Dytf=Dy_transp_3d(f,H)
M=H.vol_size(1); N=H.vol_size(2); R=H.vol_size(3);
f=reshape(f,M,N,R);
Dytf=zeros(M,N,R);

Dytf(:,N,:)=f(:,N,:);
Dytf(:,1:(N-1),:)=f(:,1:(N-1),:)-f(:,2:N,:);
Dytf=Dytf(:);