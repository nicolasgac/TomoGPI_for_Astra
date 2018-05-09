function isnr=isnr(f,fh0,fh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the Improvement Signal-to-Noise Ratio (ISNR) for
% 2D or 3D images by using iterative approaching methods.
% Input:
%       f: original image
%       fh0: initialized image
%       fh: reconstructed image
% Author: Li Wang
% January 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mse1=mean((fh(:)-f(:)).^2);
mse2=mean((fh0(:)-f(:)).^2);
isnr=10*log10(mse2/mse1);