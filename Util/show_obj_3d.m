function show_obj_3d(f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    This show_3d function shows the middle slice of phantom in four
%        slices: the 0.25, 0.4, 0.5 and 0.75 slices.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[M,N,R]=size(f);
subplot(221); imagesc(f(:,:,R/4)); colormap(gray); %colorbar;
subplot(222); imagesc(f(:,:,R/2)); colormap(gray); %colorbar;
subplot(223); imagesc(squeeze(f(:,N/2,:))); colormap(gray); %colorbar;
subplot(224); imagesc(squeeze(f(M/2,:,:))); colormap(gray); %colorbar;

return