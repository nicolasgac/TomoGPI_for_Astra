function [ssim, luminance, contrast, structure] = SSIM(F, F_hat, alpha, beta, gamma, L)
% compute the SSIM (Structural Similarity) between the true objet F and its
% estimation F_hat
% L : maximum range for pixel values (L=255 : default value, for 8-bit
% grayscale image)
% alpha, beta, gamma : positive exponents to adjust the relative importance
% of luminance, contrast and structural comparisons
%
% luminance : luminance comparison
% contrast : contrast comparison
% structure : structure comparison
%
% Reference : "Image Quality Assessment : from Error Visibility to
% Structural Similarity", Zhou Wang, Alan Conrad Bovik, Hamid Rahim Sheikh
% and Eero P. Simoncelli (2004)
%
% Author : Camille Chapdelaine
% January 2017

if nargin<3
   %default values
   L=255; 
   alpha=1;
   beta=1;
   gamma=1;
elseif nargin<7
   L=255;%default value
end

if size(F)~=size(F_hat)
   error('The two objects must be the same size') 
end

% number of pixels
N=numel(F);

% luminance comparison
mu=sum(F(:))/N;
mu_hat=sum(F_hat(:))/N;
K1=0.03;
C1=(K1*L)^2;% small constant to avoid division by zero

luminance=(2*mu*mu_hat+C1)/(mu^2+mu_hat^2+C1);

% contrast comparison
sigma2=sum((F(:)-mu).^2)/(N-1);
sigma=sqrt(sigma2);

sigma2_hat=sum((F_hat(:)-mu_hat).^2)/(N-1);
sigma_hat=sqrt(sigma2_hat);

K2=0.01;
C2=(K2*L)^2;% small constant to avoid division by zero

contrast=(2*sigma*sigma_hat+C2)/(sigma2+sigma2_hat+C2);

% structure comparison
correlation=sum((F(:)-mu).*(F_hat(:)-mu_hat))/(N-1);

C3=C2/2;%small constant to avoid division by zero

structure=(correlation+C3)/(sigma*sigma_hat+C3);

% SSIM
ssim=(luminance^alpha)*(contrast^beta)*(structure^gamma);


end

