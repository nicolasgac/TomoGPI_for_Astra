function [MeanTime,fh,ve,vxi,vz,v_erf_l2,v_erg_l2,v_erz_l2,v_erf_l1,v_erg_l1,v_erz_l1] = ...
                                     TomoGPI_HHBM_St_3d(H,g,itermax,fh0,f,l,snr,niter_gradient,az0,bz0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function uses the HHBM method where all the variances ve, vxi and vz
% are vectors. This function is used for the 3D object reconstructions
% where ASTRA toolbox is used.
% Author: Dumitru Mircea
% August 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TimeVector=[];
fh=fh0;
M=H.vol_size(1); N=H.vol_size(2); R=H.vol_size(3); MNR=M*N*R;
NS=H.proj_size(1); NP=H.proj_size(2); ND=H.proj_size(3); NSNDNP=NS*ND*NP;

% haar
z=haar3_GPU(f,M,N,R,l);

normf=sum(f(:).^2);
normg=sum(g(:).^2);
normz=sum(z(:).^2);

% haar
zh=haar3_GPU(fh,M,N,R,l);
zh=reshape(zh,M,N,R);
gh=H*fh(:);gh=reshape(gh,H.proj_size);
deltag=g(:)-gh(:);
deltaf=f(:)-fh(:);
deltaz=abs(z(:)-zh(:));
dg_l1=sum(abs(deltag))/sum(abs(g(:)));df_l1=sum(abs(deltaf))/sum(abs(f(:)));dz_l1=sum(abs(deltaz))/sum(abs(z(:)));
dg_l2=sum(deltag.^2)/normg;df_l2=sum(deltaf.^2)/normf;dz_l2=sum(deltaz.^2)/normz;
%disp('|f-fh|^2/|f|^2, |g-gh|^2/|g|^2, |z-zh|^2/|z|^2')
%disp([df_l2,dg_l2,dz_l2])
%disp('|f-fh|/|f|, |g-gh|/|g|, |z-zh|/|z|')
%disp([df_l1,dg_l1,dz_l1])

v_erf_l1=zeros(itermax+1,1);
v_erg_l1=zeros(itermax+1,1);
v_erz_l1=zeros(itermax+1,1);
v_erf_l2=zeros(itermax+1,1);
v_erg_l2=zeros(itermax+1,1);
v_erz_l2=zeros(itermax+1,1);
v_erf_l1(1)=df_l1;v_erg_l1(1)=dg_l1;v_erz_l1(1)=dz_l1;
v_erf_l2(1)=df_l2;v_erg_l2(1)=dg_l2;v_erz_l2(1)=dz_l2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ae0=2.01*ones(NSNDNP,1);
switch snr
    case 10
        be0=10*ones(NSNDNP,1); %% Very Important !!!!!!
    case 20
        be0=7*ones(NSNDNP,1); %% Very Important !!!!!!
    case 30
        be0=5*ones(NSNDNP,1); %% Very Important !!!!!!
    case 40
        be0=2*ones(NSNDNP,1); %% Very Important !!!!!!
    case 1000
        be0=2*ones(NSNDNP,1); %% Very Important !!!!!!
end

axi0=2.01*ones(MNR,1); bxi0=0.1*ones(MNR,1);

%az0=0.001;
%bz0=0.001;

%niter_gradient=20;
niter_total=itermax*niter_gradient;

critf.J=zeros(niter_total,1);
critf.J_MC=zeros(niter_total,1);
critf.J_reg=zeros(niter_total,1);
critf.alpha=zeros(niter_total,1);
critf.num_iter=0;

critz.J=zeros(niter_total,1);
critz.J_MC=zeros(niter_total,1);
critz.J_reg=zeros(niter_total,1);
critz.alpha=zeros(niter_total,1);
critz.num_iter=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Starting iterations for St prior')
disp('...')

h = waitbar(0,'Please wait...');

s = clock;
for iter=1:1:itermax
    tic
    num_iter_global=iter;
    % vz
    vz=(bz0+0.5*zh(:).^2)./(az0+1.5);
    
    % ve
    deltag=g(:)-H*fh(:);
    ve=(be0+0.5*deltag.^2)./(ae0+1.5);
    
    % vxi
    Dzh=ihaar3_GPU(zh,M,N,R,l);
    er_f=fh(:)-Dzh(:);
    vxi=(bxi0+0.5*er_f.^2)./(axi0+1.5);
    
    % f
    [fh,critf]=TomoGPI_F_Gradiant(H,fh,g,zh,num_iter_global,niter_gradient,critf,ve,vxi,M,N,R,l);
    fh(fh<0)=0;
    
    % z
    [zh,critz]=TomoGPI_Z_Gradiant(H,g,zh,fh,num_iter_global,niter_gradient,critz,vxi,ve,vz,M,N,R,l,z);
    TimeVector=[TimeVector,toc];
    %%%
    fh=reshape(fh,H.vol_size);
    zh=reshape(zh,H.vol_size);
    %fh(fh<0)=0;
    gh=H*fh(:);gh=reshape(gh,H.proj_size);
    deltag=g(:)-gh(:);
    deltaf=f(:)-fh(:);
    deltaz=z(:)-zh(:);
    dg_l1=sum(abs(deltag))/sum(abs(g(:)));df_l1=sum(abs(deltaf))/sum(abs(f(:)));dz_l1=sum(abs(deltaz))/sum(abs(z(:)));
    dg_l2=sum(deltag.^2)/normg;df_l2=sum(deltaf.^2)/normf;dz_l2=sum(abs(deltaz).^2)/normz;
     %disp('|f-fh|^2/|f|^2, |g-gh|^2/|g|^2, |z-zh|^2/|z|^2')
     %disp([df_l2,dg_l2,dz_l2])
     %disp('|f-fh|/|f|, |g-gh|/|g|, |z-zh|/|z|')
     %disp([df_l1,dg_l1,dz_l1])
    
%     figure(4);clf,show_obj_3d(fh);suptitle('f reconstructed');
%     figure(5);clf,show_proj_3d(gh);suptitle('g reconstructed');
%     figure(6);clf,imagesc(zh(:,:,R/(2^(l+1))));colormap(gray);title('z reconstructed');axis('square'),axis off;drawnow
    
    v_erf_l1(iter+1)=df_l1;
    v_erg_l1(iter+1)=dg_l1;
    v_erz_l1(iter+1)=dz_l1;
    v_erf_l2(iter+1)=df_l2;
    v_erg_l2(iter+1)=dg_l2;
    v_erz_l2(iter+1)=dz_l2;
    
    %subplot(221); imagesc(f(:,:,R/4)); colormap(gray); %colorbar;
    %subplot(222); imagesc(f(:,:,R/2)); colormap(gray); %colorbar;
    %subplot(223); imagesc(squeeze(f(:,N/2,:))); colormap(gray); %colorbar;
    %subplot(224); imagesc(squeeze(f(M/2,:,:))); colormap(gray); %colorbar;

    figure(100)    
    
        subplot(2,3,1);
        plot(v_erf_l2(1:iter+1));
        title('relative error l2 of f');
    
        subplot(2,3,4); 
        plot(v_erf_l1(1:iter+1));
        title('relative error l1 of f');
    
        subplot(2,3,2); 
        plot(v_erg_l2(1:iter+1));
        title('relative error l2 of g');
    
        subplot(2,3,5);
        plot(v_erg_l1(1:iter+1));
        title('relative error l1 of g');
    
        subplot(2,3,3); 
        plot(v_erz_l2(1:iter+1));   
        title('relative error l2 of z');
    
        subplot(2,3,6); 
        plot(v_erz_l1(1:iter+1));   
        title('relative error l1 of z');    
        
        
        if iter ==1

      is = etime(clock,s);

      esttime = is * itermax;

     end

     h = waitbar(iter/itermax,h,['remaining time =',num2str(esttime-etime(clock,s),'%4.1f'),'sec']);

end
close(h);

% figure(4);clf,plot(v_erf_l1(1:(iter+1)));title('relative error l1 of f');drawnow
% figure(5);clf,plot(v_erg_l1(1:(iter+1)));title('relative error l1 of g');drawnow
% figure(7);clf,plot(v_erf_l2(1:(iter+1)));title('relative error l2 of f');drawnow
% figure(8);clf,plot(v_erg_l2(1:(iter+1)));title('relative error l2 of g');drawnow
%     figure(9);clf,plot(v_erz_l2(1:iter+1));title('relative error of z');drawnow


fh=reshape(fh,H.vol_size);
zh=reshape(zh,H.vol_size);
gh=H*fh(:);gh=reshape(gh,H.proj_size);

 deltag=abs(g(:)-gh(:));
 deltaf=abs(f(:)-fh(:));
 deltaz=abs(z(:)-zh(:));
 dg_l1=sum(abs(deltag))/sum(abs(g(:)));df_l1=sum(abs(deltaf))/sum(abs(f(:)));dz_l1=sum(abs(deltaz))/sum(abs(z(:)));
 dg_l2=sum(deltag.^2)/normg;df_l2=sum(deltaf.^2)/normf;dz_l2=sum(deltaz.^2)/normz;
 %disp('|f-fh|^2/|f|^2, |g-gh|^2/|g|^2, |z-zh|^2/|z|^2')
 %disp([df_l2,dg_l2,dz_l2])
 %disp('|f-fh|/|f|, |g-gh|/|g|, |z-zh|/|z|')
 %disp([df_l1,dg_l1,dz_l1])

% figure(4);clf,show_obj_3d(fh);suptitle('f reconstructed');
% figure(5);clf,show_proj_3d(gh);suptitle('g reconstructed');
% figure(6);clf,imagesc(zh(:,:,R/(2^(l+1))));colormap(gray);title('z reconstructed');axis('square'),axis off;drawnow
MeanTime=mean(TimeVector);
disp('Finished iterations for St prior');
disp(['Elepsed time: ',num2str(MeanTime)]);
end
