function v=ihaar3_GPU(u,M,N,R,l)

%% ihaar3_GPU inverts the Haar transform of an 3D tensor by using GPU.
% Parameters:
%  Input, integer M, N, R the dimensions of the object.
%  Input, real U(M,N,R), the object to be transformed.
%  Output, real V(M,N,R), the transformed tensor.
%
s=sqrt(2)/2;

u=reshape(u,M,N,R);
v=u;

km=M/(2^(l-1));
kn=N/(2^(l-1));
kr=R/(2^(l-1));

while(km<=M && kn<=N && kr<=R)
    w=v;
    v(1:km,1:kn,1:2:kr)=s*(w(1:km,1:kn,1:(kr/2))+w(1:km,1:kn,(kr/2+1):kr));
    v(1:km,1:kn,2:2:kr)=s*(w(1:km,1:kn,1:(kr/2))-w(1:km,1:kn,(kr/2+1):kr));
    w=v;
    v(1:km,1:2:kn,1:kr)=s*(w(1:km,1:(kn/2),1:kr)+w(1:km,(kn/2+1):kn,1:kr));
    v(1:km,2:2:kn,1:kr)=s*(w(1:km,1:(kn/2),1:kr)-w(1:km,(kn/2+1):kn,1:kr));
    w=v;
    v(1:2:km,1:kn,1:kr)=s*(w(1:(km/2),1:kn,1:kr)+w((km/2+1):km,1:kn,1:kr));
    v(2:2:km,1:kn,1:kr)=s*(w(1:(km/2),1:kn,1:kr)-w((km/2+1):km,1:kn,1:kr));
    km=2*km;
    kn=2*kn;
    kr=2*kr;
end