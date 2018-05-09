function [fh,V_erf_l2,V_erg_l2,V_erf_l1,V_erg_l1,V_isnr,V_psnr,V_ssim,V_time]=TomoGPI_TV_3d(H,g,itermax,fh0,f,lambda,threshold,niter_gradient)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Total Variation method 
%
% Author: Mircea Dumitru and Li Wang
% October 2018, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fh=fh0;
Dxfh=Dx_direct_3d(fh,H);
Dyfh=Dy_direct_3d(fh,H);
Dzfh=Dz_direct_3d(fh,H);
dxh=shrink(Dxfh(:),threshold);
dyh=shrink(Dyfh(:),threshold);
dzh=shrink(Dzfh(:),threshold);
bxh=zeros(size(fh));bxh=bxh(:);
byh=zeros(size(fh));byh=byh(:);
bzh=zeros(size(fh));bzh=bzh(:);
gh=H*fh(:);
gh=reshape(gh,H.proj_size);
% figure(4);imagesc(fh);colormap(gray);axis('square');title('The reconstructed fh');drawnow
% figure(5);imagesc(gh);colormap(gray);axis('square');title('The reconstructed gh');drawnow

normf=sum(abs(f(:)).^2);normg=sum(g(:).^2);
V_erf_l1=zeros(itermax+1,1);V_erg_l1=zeros(itermax+1,1);
V_erf_l2=zeros(itermax+1,1);V_erg_l2=zeros(itermax+1,1);
V_isnr=zeros(itermax+1,1); V_psnr=zeros(itermax+1,1); V_ssim=zeros(itermax+1,1);
V_time=zeros(itermax);

deltf=abs(f(:)-fh(:));deltg=abs(g(:)-gh(:));
df_l1=sum(deltf(:))/sum(abs(f(:)));dg_l1=sum(deltg(:))/sum(abs(g(:)));
df_l2=sum(deltf(:).^2)/normf;dg_l2=sum(deltg(:).^2)/normg;
V_erf_l1(1)=df_l1;V_erg_l1(1)=dg_l1;
V_erf_l2(1)=df_l2;V_erg_l2(1)=dg_l2;
V_isnr(1)=isnr(f(:),fh0,fh); V_psnr(1)=psnr(f(:),fh); V_ssim(1)=SSIM(f(:),fh,1,1,1,1);

alpha=0.00001;
%alpha=1;
disp('Starting iterations for TV')
disp('...')

h = waitbar(0,'Please wait...');

s = clock;

%for iter=1:1:itermax    
   
    for it=1:1:niter_gradient
        tic
        deltg=g(:)-H*fh(:);
        dJ_MC=-2*H'*deltg;
        Dxfh=Dx_direct_3d(fh,H);
        Dyfh=Dy_direct_3d(fh,H);
        Dzfh=Dz_direct_3d(fh,H);
        deltDfh_x=dxh-bxh-Dxfh;
        deltDfh_y=dyh-byh-Dyfh;
        deltDfh_z=dzh-bzh-Dzfh;
        DtdeltDfh_x=Dx_transp_3d(deltDfh_x,H);
        DtdeltDfh_y=Dy_transp_3d(deltDfh_y,H);
        DtdeltDfh_z=Dz_transp_3d(deltDfh_z,H);
        dJ_reg=-lambda*(DtdeltDfh_x+DtdeltDfh_y+DtdeltDfh_z);
        dJ=dJ_MC+dJ_reg;
        dJ=reshape(dJ,H.vol_size);
        direction=-dJ;
        
%         alpha=Backtracking_TV(H,g,f,bxh,byh,bzh,dxh,dyh,dzh,dJ,alpha,lambda);
%         num_alpha=sum(dJ.^2);
%         HdJ=direct(H,dJ);
%         DdJ=D_direct(dJ);
%         denum_alpha=sum(HdJ.^2)+lambda*sum(DdJ.^2);
%         alpha=0.5*num_alpha/denum_alpha;
        
        fh=fh+alpha*direction(:);
   % end
    
    Dxfh=Dx_direct_3d(fh,H);
    Dyfh=Dy_direct_3d(fh,H);
    Dzfh=Dz_direct_3d(fh,H);
    dxh=shrink((Dxfh+bxh),threshold);
    dyh=shrink((Dyfh+byh),threshold);
    dzh=shrink((Dzfh+bzh),threshold);
    bxh=bxh+Dxfh-dxh;
    byh=byh+Dyfh-dyh;
    bzh=bzh+Dzfh-dzh;
    
    gh=H*fh(:);
    gh=reshape(gh,H.proj_size);
    
    deltf=abs(f(:)-fh(:));deltg=abs(g(:)-gh(:));
    df_l1=sum(deltf(:))/sum(abs(f(:)));dg_l1=sum(deltg(:))/sum(abs(g(:)));
    df_l2=sum(deltf(:).^2)/normf;dg_l2=sum(deltg(:).^2)/normg;
    V_erf_l1(it+1)=df_l1;V_erg_l1(it+1)=dg_l1;
    V_erf_l2(it+1)=df_l2;V_erg_l2(it+1)=dg_l2;
    V_isnr(it+1)=isnr(f(:),fh0,fh); V_psnr(it+1)=psnr(f(:),fh); V_ssim(it+1)=SSIM(f(:),fh,1,1,1,1);
    
     %disp('[df_l1,dg_l1]');disp([df_l1,dg_l1]);
     %disp('[df_l2,dg_l2]');disp([df_l2,dg_l2]);
%     figure(4);imagesc(fh);colormap(gray);axis('square');title('The reconstructed fh');drawnow
%     figure(5);imagesc(gh);colormap(gray);axis('square');title('The reconstructed gh');drawnow
%     figure(4);plot(V_erf_l1(1:(iter+1)));title('\delta_f l1');drawnow
%     figure(5);plot(V_erg_l1(1:(iter+1)));title('\delta_g l1');drawnow
%     figure(7);plot(V_erf_l2(1:(iter+1)));title('\delta_f l2');drawnow
%     figure(8);plot(V_erg_l2(1:(iter+1)));title('\delta_g l2');drawnow


    figure(99);    
    
        subplot(2,2,1);
        plot(V_erf_l2(1:it+1));
        title('relative error l2 of f');
    
        subplot(2,2,3); 
        plot(V_erf_l1(1:it+1));
        title('relative error l1 of f');
    
        subplot(2,2,2); 
        plot(V_erg_l2(1:it+1));
        title('relative error l2 of g');
    
        subplot(2,2,4);
        plot(V_erg_l1(1:it+1));
        title('relative error l1 of g');

V_time(it)=toc;
if it ==1

      is = etime(clock,s);

      esttime = is * niter_gradient;

     end

     h = waitbar(it/niter_gradient,h,['remaining time =',num2str(esttime-etime(clock,s),'%4.1f'),'sec']);

end
disp('Finished iterations for St prior');
fh=reshape(fh,H.vol_size);

close(h);














