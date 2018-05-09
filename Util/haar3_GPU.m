function Haar3 = haar3_GPU (u,M,N,R,l)

%*****************************************************************************
%
%  Haar3_GPU computes the Haar transform of an 3D tensor by using GPU.
%
%  Author: Li Wang
%  2016/04/14
%
%  Input, real u(M,N,R), the 3D tensor to be transformed.
%  Output, real Haar(M,N,R), the transformed 3D tensor.
%
u=reshape(u,M,N,R);
v=u;
km=M;
kn=N;
kr=R;

min_kM=M/(2^(l-1));
min_kN=N/(2^(l-1));
min_kR=R/(2^(l-1));

while(km>=min_kM && kn>=min_kN && kr>=min_kR)
    [v(1:km,1:kn,1:kr)]=cat(1, (v(1:2:km,1:kn,1:kr)+v(2:2:km,1:kn,1:kr))/sqrt(2), (v(1:2:km,1:kn,1:kr)-v(2:2:km,1:kn,1:kr))/sqrt(2));
    [v(1:km,1:kn,1:kr)]=cat(2, (v(1:km,1:2:kn,1:kr)+v(1:km,2:2:kn,1:kr))/sqrt(2), (v(1:km,1:2:kn,1:kr)-v(1:km,2:2:kn,1:kr))/sqrt(2));
    [v(1:km,1:kn,1:kr)]=cat(3, (v(1:km,1:kn,1:2:kr)+v(1:km,1:kn,2:2:kr))/sqrt(2), (v(1:km,1:kn,1:2:kr)-v(1:km,1:kn,2:2:kr))/sqrt(2));
    km=km/2;
    kn=kn/2;
    kr=kr/2;
end

Haar3=v;