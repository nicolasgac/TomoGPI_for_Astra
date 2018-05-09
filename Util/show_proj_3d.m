function show_proj_3d(g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    This show_3d function shows the middle slice of phantom in four
%        slices: the 0.25, 0.4, 0.5 and 0.75 slices.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[NS,ND,NPHI]=size(g);
subplot(121); imagesc(g(:,:,NPHI/2)); colormap(gray); axis('square'); axis off;
subplot(122); imagesc(squeeze(g(:,ND/2,:))); colormap(gray); axis('square'); axis off;

return