%---------------------------------------------------------------------------------------------------%
%          3D object reconstruction using ASTRA toolbox and Spot operator                           %
%---------------------------------------------------------------------------------------------------%
% Menu for demonstrative purposes;                                                                  %
% Testing / showing bayesian algorithms in CT for ASTRA library integration                         %
%---------------------------------------------------------------------------------------------------%
% FUNCTIONS USED in this program :                                                                  %
%                                                                                                   %
%    ASTRA functions :                                                                              %
%        astra_create_vol_geom                                                                      %
%        astra_create_proj_geom                                                                     %
%        opTomo                                                                                     %
%        show_proj_3d                                                                               %
%        show_obj_3d                                                                                %
%                                                                                                   %
%    MATLAB build in functions :                                                                    %
%        linspace2                                                                                  %
%                                                                                                   %
%    GPI functions :                                                                                %
%        Haar_JMAP_St_3d  (Bayesian approach, Student-t prior, 3d dedicated)                        %
%        Haar_JMAP_NIG_3d (Bayesian approach, Normal Inverse Gaussian prior, 3d dedicated)          %
%        Haar_JMAP_VG_3d  (Bayesian approach, Variance Gamma prior, 3d dedicated)                   %
%                ihaar3_GPU                                                                         %
%---------------------------------------------------------------------------------------------------%
% nproj:  number of projections                                                                     %
% snr:    signal to noise ratio, used to define ve: ve=vg/(10^(snr/10))                             %
% phantom3d defines the 3D Shepp Logan object                                                       %
%---------------------------------------------------------------------------------------------------%
% Results are stored in .mat files for each launch of reconstruction (FBP, TV,
% St, NIG or VG) stamped with the date of launch
% to load them use "test1_TV=load(./TV/test_date_parameters.mat)"
% them access it through "test1_TV.variable" for instance test1_TV.fhTV  
%---------------------------------------------------------------------------------------------------%
% Authors            :  Mircea Dumitru (uses code previusly developed by Contributors)  &  Nicolas Gac    %
% Supervisers       :  Nicolas Gac  &  Ali Mohammad-Djafari                                         %
% Contributors      :  Li Wang                                                                      %
% April 2018,  @ GPI,  L2S, CentraleSupelec <- last modifications                                   %
%---------------------------------------------------------------------------------------------------%
%---------------------------------------------------------------------------------------------------%
clear all;          close all;          clc;                                                        %
%---------------------------------------------------------------------------------------------------%
% default values;                                                                                   %
% if the user selects to run the algorithm directly, default valus are used;                        %
%---------------------------------------------------------------------------------------------------%
save_path       =  '.';
path_menu=evalc('which Menu_TomoGPI_ASTRA');
[pathstr,name,ext] = fileparts(path_menu);
addpath([pathstr , '/Util']);
%---------------------------------------------------------------------------------------------------%
fn_FBP          =  'FBP';                                                                           %
save_path_FBP   =  [save_path, '/', fn_FBP];                                                        %
mkdir(save_path,fn_FBP); 
counter_FBP=1;%
%---------------------------------------------------------------------------------------------------%
fn_TV           =  'TV';                                                                            %
save_path_TV    =  [save_path, '/', fn_TV];                                                         %
mkdir(save_path,fn_TV);   
counter_TV=1;%
%---------------------------------------------------------------------------------------------------%
fn_St           =  'St';                                                                            %
save_path_St    =  [save_path, '/', fn_St];                                                         %
mkdir(save_path,fn_St);    
counter_St=1;%
%---------------------------------------------------------------------------------------------------%
fn_NIG          =  'NIG';                                                                           %
save_path_NIG   =  [save_path, '/', fn_NIG];                                                        %
mkdir(save_path,fn_NIG);  
counter_NIG=1;%
%---------------------------------------------------------------------------------------------------%
fn_VG           =  'VG';                                                                            %
save_path_VG    =  [save_path, '/', fn_VG];                                                         %
mkdir(save_path,fn_VG);   
counter_VG=1;%
%---------------------------------------------------------------------------------------------------%
% Default values if the user is actually skiping the buttons where it sets                          %
% data and goes directly to the method.                                                             %
%---------------------------------------------------------------------------------------------------%
N               =  64;                                                                              %
nproj           =  64;                                                                              %
snr             =  30;                                                                              %
%--------------------                                                                               %
f               =  phantom3d(N); % function phantom3d creats a 3D phantom, size N^3;                %
vol_geom        =  astra_create_vol_geom(N, N, N);                                                  %
proj_geom       =  astra_create_proj_geom('parallel3d', 1, 1, N, N, linspace2(0,pi,nproj));         %
H_base          =  opTomo('cuda', proj_geom, vol_geom);                                             %
g               =  H_base * f(:);                                                                   %
g               =  awgn(g, snr, 'measured');                                                        %
g               =  reshape(g, H_base.proj_size);                                                    %
%--------------------                                                                               %
fh0_base        =  H_base' * g(:);                                                                  %
gh0_base        =  H_base * fh0_base;                                                               %
fact_base       =  norm(g(:)) / norm(gh0_base(:));                                                  %
fh0_base        =  fact_base * fh0_base;                                                            %
%---------------------------------------------------------------------------------------------------%
% Default values for the parameters corresponding to different methods.                             %
%---------------------------------------------------------------------------------------------------%
itermax         =  20;                 % Bayesian                                                   %
niter_gradient  =  20;                 % Bayesian                                                   %
%---------------------------------------------------------------------------------------------------%
azSt0           =  2.1;                % Bayesian St                                                %
bzSt0           =  2.1;                % Bayesian St                                                %
%---------------------------------------------------------------------------------------------------%
azVG0           =  2.1;                % Bayesian VG                                                %
bzVG0           =  1/2.1;              % Bayesian VG                                                %
%---------------------------------------------------------------------------------------------------%
azNIG0          =  2.1;                % Bayesian NIG                                               %
bzNIG0          =  2.1;                % Bayesian NIG                                               %
%---------------------------------------------------------------------------------------------------%
l               =  5;                  % Bayesian                                                   %
%---------------------------------------------------------------------------------------------------%
lambda          =  1;                  % TV                                                         %
threshold       =  1;                  % TV                                                         %
err_min         =  1;                  % TV                                                         %
%---------------------------------------------------------------------------------------------------%
sur_ech         =  16;                 % FBK                                                        %
fc_FDK          =  1;                  % FBK                                                        %
delta_un        =  1;                  % FBK                                                        %
%---------------------------------------------------------------------------------------------------%
h_2d            =  [0 -1 0; -1 4 -1; 0 -1 0];                                                       %
%---------------------------------------------------------------------------------------------------%
suffix          =  strcat(int2str(nproj), 'proj_', int2str(snr), 'dB');                             %
%---------------------------------------------------------------------------------------------------%
while 1                                                                                             %
choice  =  menu('Menu', '1. Chose data', '2. Show Data', '3. Method', '4. Quit');                   %
switch choice                                                                                       %
%---------------------------------------------------------------------------------------------------%
%        Data button                                                                                %
%---------------------------------------------------------------------------------------------------%
case 1                                                                                              %
    while 1                                                                                         %
    choiceDataProj       =  menu('Data & Projection & SNR', 'SNR (default 30)',...                  %
                                                            'Data dimension (default 64)',...       %
                                                            'Projection (default 64)','Back');      %
    switch choiceDataProj                                                                           %
    %-----------------------------------------------------------------------------------------------%
    case 1                                                                                          %
        while 1                                                                                     %
        choiceSNR        =  menu('SNR', '40', '30', '20', '10');                                    %
        switch choiceSNR                                                                            %
        case 1                                                                                      %
            snr          =  40;                                                                     %
            str=['Setting the noise level to SNR = ', num2str(snr)]; disp(str)                      %
            str=['....']; disp(str)                                                                 %
            suffix       =  strcat(int2str(nproj),'proj_',int2str(snr),'dB');                       %
            %--------------------                                                                   %
            str=['The noise level is now set to SNR = ', num2str(snr)]; disp(str)                   %
            str=['.....[ The data dimension was set to N = ', num2str(N),' ]'];                     %
            disp(str)                                                                               %
            str=['.....[ The projection number was set to nproj = ', num2str(nproj),' ]'];          %
            disp(str)                                                                               %
            str=['   ']; disp(str)                                                                  %
            break;                                                                                  %
        case 2                                                                                      %
            snr          =  30;                                                                     %
            str=['Setting the noise level to SNR = ', num2str(snr)]; disp(str)                      %
            str=['....']; disp(str)                                                                 %
            suffix       =  strcat(int2str(nproj),'proj_',int2str(snr),'dB');                       %
            %--------------------                                                                   %
            str=['The noise level is now set to SNR = ', num2str(snr)]; disp(str)                   %
            str=['.....[ The data dimension was set to N = ', num2str(N),' ]'];                     %
            disp(str)                                                                               %
            str=['.....[ The projection number was set to nproj = ', num2str(nproj),' ]'];          %
            disp(str)                                                                               %
            str=['   ']; disp(str)                                                                  %
            break;                                                                                  %
        case 3                                                                                      %
            snr          =  20;                                                                     %
            str=['Setting the noise level to SNR = ', num2str(snr)]; disp(str)                      %
            str=['....']; disp(str)                                                                 %
            suffix       =  strcat(int2str(nproj),'proj_',int2str(snr),'dB');                       %
            %--------------------                                                                   %
            str=['The noise level is now set to SNR = ', num2str(snr)]; disp(str)                   %
            str=['.....[ The data dimension was set to N = ', num2str(N),' ]'];                     %
            disp(str)                                                                               %
            str=['.....[ The projection number was set to nproj = ', num2str(nproj),' ]'];          %
            disp(str)                                                                               %
            str=['   ']; disp(str)                                                                  %
            break;                                                                                  %
        case 4                                                                                      %
            snr          =  10;                                                                     %
            str=['Setting the noise level to SNR = ', num2str(snr)]; disp(str)                      %
            str=['....']; disp(str)                                                                 %
            suffix       =  strcat(int2str(nproj),'proj_',int2str(snr),'dB');                       %
            %--------------------                                                                   %
            str=['The noise level is now set to SNR = ', num2str(snr)]; disp(str)                   %
            str=['.....[ The data dimension was set to N = ', num2str(N),' ]'];                     %
            disp(str)                                                                               %
            str=['.....[ The projection number was set to nproj = ', num2str(nproj),' ]'];          %
            disp(str)                                                                               %
            str=['   ']; disp(str)                                                                  %
            break;                                                                                  %
        case 5                                                                                      %
            break;                                                                                  %
        end                                                                                         %
        end                                                                                         %
    %-----------------------------------------------------------------------------------------------%
    case 2                                                                                          %
        while 1                                                                                     %
        choiceData       =  menu('Data ShepLog 3D dimension', '1024', '256', '128', '64', 'Back');  %
        switch choiceData                                                                           %
        case 1                                                                                      %
            N            =  1024;                                                                   %
            str=['Setting data dimension to N = ', num2str(N), ' [this might take a few seconds]']; %
            disp(str)                                                                               %
            str=['....']; disp(str)                                                                 %
            f            =  phantom3d(N); % function phantom3d creats a 3D phantom, size N^3;       %
            vol_geom     =  astra_create_vol_geom(N, N, N);                                         %
            proj_geom    =  astra_create_proj_geom('parallel3d', 1, 1, N, N, linspace2(0,pi,nproj));%
            H_base       =  opTomo('cuda', proj_geom, vol_geom);                                    %
            g            =  H_base * f(:);                                                          %
            g            =  awgn(g,snr,'measured');                                                 %
            g            =  reshape(g, H_base.proj_size);                                           %
            str=['Data dimension is now set to N = ', num2str(N)]; disp(str)                        %
            str=['.....[ The number of projections is (by default) nproj = ', num2str(nproj),' ]']; %
            disp(str)                                                                               %
            str=['.....[ The level of noise is (by default) SNR = ', num2str(snr),' ]'];            %
            disp(str)                                                                               %
            str=['   ']; disp(str)                                                                  %
            break;                                                                                  %
        case 2                                                                                      %
            N            =  256;                                                                    %
            str=['Setting data dimension to N = ', num2str(N)]; disp(str)                           %
            str=['....']; disp(str)                                                                 %
            f            =  phantom3d(N); % function phantom3d creats a 3D phantom, size N^3;       %
            vol_geom     =  astra_create_vol_geom(N, N, N);                                         %
            proj_geom    =  astra_create_proj_geom('parallel3d', 1, 1, N, N, linspace2(0,pi,nproj));%
            H_base       =  opTomo('cuda', proj_geom, vol_geom);                                    %
            g            =  H_base * f(:);                                                          %
            g            =  awgn(g,snr,'measured');                                                 %
            g            =  reshape(g, H_base.proj_size);                                           %
            str=['Data dimension is now set to N = ', num2str(N)]; disp(str)                        %
            str=['.....[ The number of projections is (by default) nproj = ', num2str(nproj),' ]']; %
            disp(str)                                                                               %
            str=['.....[ The level of noise is (by default) SNR = ', num2str(snr),' ]'];            %
            disp(str)                                                                               %
            str=['   ']; disp(str)                                                                  %
            break;                                                                                  %
        case 3                                                                                      %
            N            =  128;                                                                    %
            str=['Setting data dimension to N = ', num2str(N)]; disp(str)                           %
            str=['....']; disp(str)                                                                 %
            f            =  phantom3d(N); % function phantom3d creats a 3D phantom, size N^3;       %
            vol_geom     =  astra_create_vol_geom(N, N, N);                                         %
            proj_geom    =  astra_create_proj_geom('parallel3d', 1, 1, N, N, linspace2(0,pi,nproj));%
            H_base       =  opTomo('cuda', proj_geom, vol_geom);                                    %
            g            =  H_base * f(:);                                                          %
            g            =  awgn(g,snr,'measured');                                                 %
            g            =  reshape(g, H_base.proj_size);                                           %
            str=['Data dimension is now set to N = ', num2str(N)]; disp(str)                        %
            str=['.....[ The number of projections is (by default) nproj = ', num2str(nproj),' ]']; %
            disp(str)                                                                               %
            str=['.....[ The level of noise is (by default) SNR = ', num2str(snr),' ]'];            %
            disp(str)                                                                               %
            str=['   ']; disp(str)                                                                  %
            break;                                                                                  %
        case 4                                                                                      %
            N            =  64;                                                                     %
            str=['Setting data dimension to N = ', num2str(N)]; disp(str)                           %
            str=['....']; disp(str)                                                                 %
            f            =  phantom3d(N); % function phantom3d creats a 3D phantom, size N^3;       %
            vol_geom     =  astra_create_vol_geom(N, N, N);                                         %
            proj_geom    =  astra_create_proj_geom('parallel3d', 1, 1, N, N, linspace2(0,pi,nproj));%
            H_base       =  opTomo('cuda', proj_geom, vol_geom);                                    %
            g            =  H_base * f(:);                                                          %
            g            =  awgn(g,snr,'measured');                                                 %
            g            =  reshape(g, H_base.proj_size);                                           %
            str=['Data dimension is now set to N = ', num2str(N)]; disp(str)                        %
            str=['.....[ The number of projections is (by default) nproj = ', num2str(nproj),' ]']; %
            disp(str)                                                                               %
            str=['.....[ The level of noise is (by default) SNR = ', num2str(snr),' ]'];            %
            disp(str)                                                                               %
            str=['   ']; disp(str)                                                                  %
            break;                                                                                  %
        case 5                                                                                      %
            break;                                                                                  %
        end                                                                                         %
        end                                                                                         %
    %-----------------------------------------------------------------------------------------------%
    case 3                                                                                          %
        while 1                                                                                     %
        choiceNProj      =  menu('No. Proj', '1024', '256', '128', '64', '32', 'Input', 'Back');    %
        switch choiceNProj                                                                          %
        case 1                                                                                      %
            nproj        =  1024;                                                                   %
            str=['Setting the projection number to nproj = ', num2str(nproj)]; disp(str)            %
            str=['....']; disp(str)                                                                 %
            %-----                                                                                  %
            vol_geom     =  astra_create_vol_geom(N, N, N);                                         %
            proj_geom    =  astra_create_proj_geom('parallel3d', 1, 1, N, N, linspace2(0,pi,nproj));%
            H_base       =  opTomo('cuda', proj_geom, vol_geom);                                    %
            g            =  H_base * f(:);                                                          %
            g            =  awgn(g,snr,'measured');                                                 %
            g            =  reshape(g, H_base.proj_size);                                           %
            %-----                                                                                  %
            fh0_base     =  H_base' * g(:);                                                         %
            gh0_base     =  H_base * fh0_base;                                                      %
            fact_base    =  norm(g(:)) / norm(gh0_base(:));                                         %
            fh0_base     =  fact_base * fh0_base;                                                   %
            %--------------------                                                                   %
            str=['The projection number is now set to nproj = ', num2str(nproj)]; disp(str)         %
            str=['.....[ The data dimension was set to N = ', num2str(N),' ]'];                     %
            disp(str)                                                                               %
            str=['.....[ The level of noise is (by default) SNR = ', num2str(snr),' ]'];            %
            disp(str)                                                                               %
            str=['   ']; disp(str)                                                                  %
            break;                                                                                  %
        %-------------------------------------------------------------------------------------------%
        case 2                                                                                      %
            nproj        =  256;                                                                    %
            str=['Setting the projection number to nproj = ', num2str(nproj)]; disp(str)            %
            str=['....']; disp(str)                                                                 %
            %-----                                                                                  %
            vol_geom     =  astra_create_vol_geom(N, N, N);                                         %
            proj_geom    =  astra_create_proj_geom('parallel3d', 1, 1, N, N, linspace2(0,pi,nproj));%
            H_base       =  opTomo('cuda', proj_geom, vol_geom);                                    %
            g            =  H_base * f(:);                                                          %
            g            =  awgn(g,snr,'measured');                                                 %
            g            =  reshape(g, H_base.proj_size);                                           %
            %-----                                                                                  %
            fh0_base     =  H_base' * g(:);                                                         %
            gh0_base     =  H_base * fh0_base;                                                      %
            fact_base    =  norm(g(:)) / norm(gh0_base(:));                                         %
            fh0_base     =  fact_base * fh0_base;                                                   %
            %--------------------                                                                   %
            str=['The projection number is now set to nproj = ', num2str(nproj)]; disp(str)         %
            str=['.....[ The data dimension was set to N = ', num2str(N),' ]'];                     %
            disp(str)                                                                               %
            str=['.....[ The level of noise is (by default) SNR = ', num2str(snr),' ]'];            %
            disp(str)                                                                               %
            str=['   ']; disp(str)                                                                  %
            break;                                                                                  %
        %-------------------------------------------------------------------------------------------%
        case 3                                                                                      %
            nproj        =  128;                                                                    %
            str=['Setting the projection number to nproj = ', num2str(nproj)]; disp(str)            %
            str=['....']; disp(str)                                                                 %
            %-----                                                                                  %
            vol_geom     =  astra_create_vol_geom(N, N, N);                                         %
            proj_geom    =  astra_create_proj_geom('parallel3d', 1, 1, N, N, linspace2(0,pi,nproj));%
            H_base       =  opTomo('cuda', proj_geom, vol_geom);                                    %
            g            =  H_base * f(:);                                                          %
            g            =  awgn(g,snr,'measured');                                                 %
            g            =  reshape(g, H_base.proj_size);                                           %
            %-----                                                                                  %
            fh0_base     =  H_base' * g(:);                                                         %
            gh0_base     =  H_base * fh0_base;                                                      %
            fact_base    =  norm(g(:)) / norm(gh0_base(:));                                         %
            fh0_base     =  fact_base * fh0_base;                                                   %
            %--------------------                                                                   %
            str=['The projection number is now set to nproj = ', num2str(nproj)]; disp(str)         %
            str=['.....[ The data dimension was set to N = ', num2str(N),' ]'];                     %
            disp(str)                                                                               %
            str=['.....[ The level of noise is (by default) SNR = ', num2str(snr),' ]'];            %
            disp(str)                                                                               %
            str=['   ']; disp(str)                                                                  %
            break;                                                                                  %
        %-------------------------------------------------------------------------------------------%
        case 4                                                                                      %
            nproj        =  64;                                                                     %
            str=['Setting the projection number to nproj = ', num2str(nproj)]; disp(str)            %
            str=['....']; disp(str)                                                                 %
            %--------------------                                                                   %
            vol_geom     =  astra_create_vol_geom(N, N, N);                                         %
            proj_geom    =  astra_create_proj_geom('parallel3d', 1, 1, N, N, linspace2(0,pi,nproj));%
            H_base       =  opTomo('cuda', proj_geom, vol_geom);                                    %
            g            =  H_base * f(:);                                                          %
            g            =  awgn(g,snr,'measured');                                                 %
            g            =  reshape(g, H_base.proj_size);                                           %
            %--------------------                                                                   %
            fh0_base     =  H_base' * g(:);                                                         %
            gh0_base     =  H_base * fh0_base;                                                      %
            fact_base    =  norm(g(:)) / norm(gh0_base(:));                                         %
            fh0_base     =  fact_base * fh0_base;                                                   %
            %--------------------                                                                   %
            str=['The projection number is now set to nproj = ', num2str(nproj)]; disp(str)         %
            str=['.....[ The data dimension was set to N = ', num2str(N),' ]'];                     %
            disp(str)                                                                               %
            str=['.....[ The level of noise is (by default) SNR = ', num2str(snr),' ]'];            %
            disp(str)                                                                               %
            str=['   ']; disp(str)                                                                  %
            break;                                                                                  %
        %-------------------------------------------------------------------------------------------%
        case 5                                                                                      %
            nproj        =  32;                                                                     %
            %--------------------                                                                   %
            str=['Setting the projection number to nproj = ', num2str(nproj)]; disp(str)            %
            str=['....']; disp(str)                                                                 %
            %--------------------                                                                   %
            vol_geom     =  astra_create_vol_geom(N, N, N);                                         %
            proj_geom    =  astra_create_proj_geom('parallel3d', 1, 1, N, N, linspace2(0,pi,nproj));%
            H_base       =  opTomo('cuda', proj_geom, vol_geom);                                    %
            g            =  H_base*f(:);                                                            %
            g            =  awgn(g,snr,'measured');                                                 %
            g            =  reshape(g, H_base.proj_size);                                           %
            %--------------------                                                                   %
            fh0_base     =  H_base' * g(:);                                                         %
            gh0_base     =  H_base * fh0_base;                                                      %
            fact_base    =  norm(g(:)) / norm(gh0_base(:));                                         %
            fh0_base     =  fact_base * fh0_base;                                                   %
            %--------------------                                                                   %
            str=['The projection number is now set to nproj = ', num2str(nproj)]; disp(str)         %
            str=['.....[ The data dimension was set to N = ', num2str(N),' ]'];                     %
            disp(str)                                                                               %
            str=['.....[ The level of noise is (by default) SNR = ', num2str(snr),' ]'];            %
            disp(str)                                                                               %
            str=['   ']; disp(str)                                                                  %
            break;                                                                                  %
        %-------------------------------------------------------------------------------------------%
        case 6                                                                                      %
            str = '[numerical values must be inputed using the keyboad, then hit enter]';           %
            disp(str)                                                                               %
            nproj        =  input('Input nproj (default 64) :');                                    %
            %--------------------                                                                   %
            str=['Setting the projection number to nproj = ', num2str(nproj)]; disp(str)            %
            str=['....']; disp(str)                                                                 %
            %--------------------                                                                   %
            vol_geom     =  astra_create_vol_geom(N, N, N);                                         %
            proj_geom    =  astra_create_proj_geom('parallel3d', 1, 1, N, N, linspace2(0,pi,nproj));%
            H_base       =  opTomo('cuda', proj_geom, vol_geom);                                    %
            g            =  H_base*f(:);                                                            %
            g            =  awgn(g,snr,'measured');                                                 %
            g            =  reshape(g, H_base.proj_size);                                           %
            %--------------------                                                                   %
            fh0_base     =  H_base' * g(:);                                                         %
            gh0_base     =  H_base * fh0_base;                                                      %
            fact_base    =  norm(g(:)) / norm(gh0_base(:));                                         %
            fh0_base     =  fact_base * fh0_base;                                                   %
            %--------------------                                                                   %
            str=['The projection number is now set to nproj = ', num2str(nproj)]; disp(str)         %
            str=['.....[ The data dimension was set to N = ', num2str(N),' ]'];                     %
            disp(str)                                                                               %
            str=['.....[ The level of noise is (by default) SNR = ', num2str(snr),' ]'];            %
            disp(str)                                                                               %
            str=['   ']; disp(str)                                                                  %
            break;                                                                                  %
        %-------------------------------------------------------------------------------------------%
        case 7                                                                                      %
            break;                                                                                  %
        end                                                                                         %
        end                                                                                         %
    %-----------------------------------------------------------------------------------------------%
    case 4                                                                                          %
        break;                                                                                      %
    end                                                                                             %
    end                                                                                             %
%---------------------------------------------------------------------------------------------------%
%        Visualisation button                                                                       %
%---------------------------------------------------------------------------------------------------%
case 2                                                                                              %
     figure(1); show_obj_3d(f);                                                                     %
     suptitle('Shepp Logan original phantom')                                                       %
     save('/espace/tomo_gpi/ResultsMircea/Original_256_SheppLogan.mat','f');                        %
     %--------------------                                                                          %
     figure(2); show_proj_3d(g);                                                                    %
     suptitle('Projection of Shepp Logan original phantom')                                         %
     save('/espace/tomo_gpi/ResultsMircea/Proj_256_SheppLogan.mat','g');                            %
     %--------------------                                                                          %
     str=['Data corresponding to dimension N = ', num2str(N), ' and projection number nproj = ',... %
          num2str(nproj), ' is ploted in Figure 1 and Figure 2.'];                                  %
     disp(str)                                                                                      %
%---------------------------------------------------------------------------------------------------%
%        Method button                                                                              %
%---------------------------------------------------------------------------------------------------%
case 3                                                                                              %
    while 1                                                                                         %
    choiceMethod         =  menu('MenuMethod', 'No. Global Iterations (default 20)', 'No. Gradient Iterations for HHBM Methods (default 20)', 'FBP','TV',...        %
                                               'HHBM - St-t', 'HHBM - NIG', 'HHBM - VG', 'Back');   %
    switch choiceMethod                                                                             %
    %-----------------------------------------------------------------------------------------------%
    case 1                                                                                          %
        while 1                                                                                     %
        choiceIterat     =  menu('No. Global Iterations (default 20) ', '30', '40', '50', '100',...        %
                                                                 'Input', 'Back');                  %
        switch choiceIterat                                                                         %
        case 1                                                                                      %
            itermax      =  30;                                                                     %
            str=['Number of iterations is now set at itermax = ', num2str(itermax)]; disp(str)      %
            break;                                                                                  %
        case 2                                                                                      %
            itermax      =  40;                                                                     %
            str=['Number of iterations is now set at itermax = ', num2str(itermax)]; disp(str)      %
            break;                                                                                  %
        case 3                                                                                      %
            itermax      =  50;                                                                     %
            str=['Number of iterations is now set at itermax = ', num2str(itermax)]; disp(str)      %
            break;                                                                                  %
        case 4                                                                                      %
            itermax      =  100;                                                                    %
            str=['Number of iterations is now set at itermax = ', num2str(itermax)]; disp(str)      %
            break;                                                                                  %
        case 5                                                                                      %
            str = '[numerical values must be inputed using the keyboad, then hit enter]';           %
            disp(str)                                                                               %
            itermax      =  input('Input nproj (default 20) :');                                    % 
            str=['Number of iterations is now set at itermax = ', num2str(itermax)]; disp(str)      %
            break;                                                                                  %
        case 6                                                                                      %
            break;                                                                                  %
        end                                                                                         %
        end
        %-----------------------------------------------------------------------------------------------%
    case 2                                                                                          %
        while 1                                                                                     %
        choiceIterat     =  menu('No. HHBM Grad. Iterations (default 20) ', '30', '40', '50', '100',...        %
                                                                 'Input', 'Back');                  %
        switch choiceIterat                                                                         %
        case 1                                                                                      %
            niter_gradient  =        30;                                                                     %
            str=['Number of iterations is now set at niter_gradient  = = ', num2str(niter_gradient)]; disp(str)      %
            break;                                                                                  %
        case 2                                                                                      %
             niter_gradient  =       40;                                                                     %
            str=['Number of iterations is now set at niter_gradient  = = ', num2str(niter_gradient)]; disp(str)      %
            break;                                                                                  %
        case 3                                                                                      %
             niter_gradient  =      50;                                                                     %
            str=['Number of iterations is now set at niter_gradient  = = ', num2str(niter_gradient)]; disp(str)      %
            break;                                                                                  %
        case 4                                                                                      %
             niter_gradient  =      100;                                                                    %
            str=['Number of iterations is now set at  niter_gradient  == ', num2str(niter_gradient)]; disp(str)      %
            break;                                                                                  %
        case 5                                                                                      %
            str = '[numerical values must be inputed using the keyboad, then hit enter]';           %
            disp(str)                                                                               %
             niter_gradient  =      input('Input nproj (default 20) :');                                    % 
            str=['Number of iterations is now set at niter_gradient  = = ', num2str(niter_gradient)]; disp(str)      %
            break;                                                                                  %
        case 6                                                                                      %
            break;                                                                                  %
        end                                                                                         %
        end        %
    case 3                                                                                          %
        while 1                                                                                     %
        choiceHyperPar   =  menu('FBP ', 'sur_ech (default 16)', 'fc_FDK (default 1)',...           %
                                         'delta_un (default 1)', 'Launch', 'Back');                 %
        switch choiceHyperPar                                                                       %
        %-------------------------------------------------------------------------------------------%
        case 1                                                                                      %
            str = '[numerical values must be inputed using the keyboad, then hit enter]';           %
            disp(str)                                                                               %
            sur_ech      =  input('Input sur_ech (default 16):');                                   %
            break;                                                                                  %
        case 2                                                                                      %
            str = '[numerical values must be inputed using the keyboad, then hit enter]';           %
            disp(str)                                                                               %
            fc_FDK       =  input('Input fc_FDK (default 1):');                                     %
            break;                                                                                  %
        case 3                                                                                      %
            str = '[numerical values must be inputed using the keyboad, then hit enter]';           %
            disp(str)                                                                               %
            delta_un     =  input('Input delta_un (default 1):');                                   %
            break;                                                                                  %
        case 4                                                                                      %
            %---------------------------------------------------------------------------------------%
            %                              Filters Back Projection (FBP)                            %
            %---------------------------------------------------------------------------------------%
            g_filtered   =  TomoGPI_FBPFilter_3d(g,sur_ech,fc_FDK,delta_un);                        %
            %---------------------------------------------------------------------------------------%
            fhFBP        =  H_base'*g_filtered(:);                                                  %
            fhFBP        =  reshape(fhFBP,H_base.vol_size);                                         %
            fhFBP        =  (fhFBP-min(fhFBP(:)))/(max(fhFBP(:))-min(fhFBP(:)));                    %
            %---------------------------------------------------------------------------------------%
            ghFBP        =  H_base*fhFBP(:); ghFBP=reshape(ghFBP,H_base.proj_size);                 %
            %---------------------------------------------------------------------------------------%
            figure(100+counter_FBP); show_obj_3d(fhFBP);                                                          %
            suptitle('Shepp Logan FBP reconstructed phantom')  
            RepName=['test_',datestr(now, 'dd-mmm-yyyy-HH-MM')];
            mkdir([save_path_FBP,'/',RepName]);
            save([save_path_FBP,'/',RepName, '/FBP_N', num2str(N), '_nProj', num2str(nproj), '_SNR',...         %
                  num2str(snr), '_surech', num2str(sur_ech), '_fcFDK', num2str(fc_FDK), '_deltaun',...%
                  num2str(delta_un), '.mat'],'fhFBP');                                   %
            % save('/espace/tomo_gpi/ResultsMircea/FBP_256_SheppLogan.mat','fhFBP');                %
            counter_FBP=counter_FBP+1;
            break;                                                                                  %
        case 5                                                                                      %
            break;                                                                                  %
        end                                                                                         %
        end                                                                                         %
    case 4                                                                                          %
        while 1                                                                                     %
        choiceHyperPar   =  menu('TV', 'lambda (default 1)', 'threshold (default 1)',...            %
                                       'Launch', 'Back');                                           %
        switch choiceHyperPar                                                                       %
        %-------------------------------------------------------------------------------------------%
        case 1                                                                                      %
            str = '[numerical values must be inputed using the keyboad, then hit enter]';           %
            disp(str)                                                                               %
            lambda       =  input('Input lambda (default 1) :');                                    %
        case 2                                                                                      %
            str = '[numerical values must be inputed using the keyboad, then hit enter]';           %
            disp(str)                                                                               %
            threshold    =  input('Input threshold (default 1) :');                                 %
        case 3                                                                                      %
            %---------------------------------------------------------------------------------------%
            %                                   Total Variation                                     %
            %---------------------------------------------------------------------------------------%
            str=['Total Variation reconstructruction launched '];                                   %
            disp(str)                                                                               %
            str=['for data with N = ', num2str(N), ', nproj = ', num2str(nproj),...                 %
                 ' and snr = ', num2str(snr)];                                                      %
            disp(str)                                                                               %
            strPrior     =  '_TV_';                                                                 %
            str=['and parameters niter = ', num2str(itermax), ', lambda = ', num2str(lambda),...    %
                 ', threshold = ', num2str(threshold),...                                           %
                 ' and niter_gradient = ', num2str(niter_gradient)];                                %
            disp(str)                                                                               %
            %---------------------------------------------------------------------------------------%
            [fhTV,V_erf_l2,V_erg_l2,V_erf_l1,V_erg_l1,V_isnr,V_psnr,V_ssim,V_time] = ...            %
                     TomoGPI_TV_3d(H_base,g,itermax,fh0_base,f,lambda,threshold,niter_gradient);    %
            ghTV         =  H_base*fhTV(:); ghTV  =  reshape(ghTV,H_base.proj_size);                %
            %---------------------------------------------------------------------------------------%
            figure(200+counter_TV); show_obj_3d(fhTV);                                                           %
            suptitle('Shepp Logan TV reconstructed phantom')  
             RepName=['test_',datestr(now, 'dd-mmm-yyyy-HH-MM')];
            %mkdir([save_path_TV,'/',RepName]);%
            filename=[save_path_TV,'/',RepName, '_TV_N', num2str(N), '_nProj', num2str(nproj), '_SNR',...           %
                  num2str(snr), '_nIter', num2str(itermax), '_lamd', num2str(lambda), '_th',...     %
                  num2str(threshold), '.mat'];
            save(filename);                                   %
            % save('/espace/tomo_gpi/ResultsMircea/TV_256_SheppLogan.mat','fhTV');                  %
            %---------------------------------------------------------------------------------------%
%             nomf_erf_l2  =  strcat('V_erf_l2_','phantom',int2str(N),strPrior,suffix);               %
%             nomf_erg_l2  =  strcat('V_erg_l2_','phantom',int2str(N),strPrior,suffix);               %
%             nomf_erf_l1  =  strcat('V_erf_l1_','phantom',int2str(N),strPrior,suffix);               %
%             nomf_erg_l1  =  strcat('V_erg_l1_','phantom',int2str(N),strPrior,suffix);               %
%             nomf_fh      =  strcat('fh_','phantom',int2str(N),strPrior,suffix);                     %
%             nomf_gh      =  strcat('gh_','phantom',int2str(N),strPrior,suffix);                     %
            %---------------------------------------------------------------------------------------%                   %
            save('-v7.3',filename,'fhTV','V_erf_l2','V_erg_l2','V_erf_l1','V_erg_l1','ghTV');    
%             save('-v7.3',[save_path_TV,'/',RepName,'/',nomf_erf_l2,'.mat'],'V_erf_l2');                         %
%             save('-v7.3',[save_path_TV,'/',RepName,'/',nomf_erg_l2,'.mat'],'V_erg_l2');                         %
%             save('-v7.3',[save_path_TV,'/',RepName,'/',nomf_erf_l1,'.mat'],'V_erf_l1');                         %
%             save('-v7.3',[save_path_TV,'/',RepName,'/',nomf_erg_l1,'.mat'],'V_erg_l1');                         %
%             save('-v7.3',[save_path_TV,'/',RepName,'/',nomf_fh,'.mat'],'fhTV');                                 %
%             save('-v7.3',[save_path_TV,'/',RepName,'/',nomf_gh,'.mat'],'ghTV');                                 %
            %---------------------------------------------------------------------------------------%  
            counter_TV=counter_TV+1;
            break;                                                                                  %
        case 4                                                                                      %
            break;                                                                                  %
        end                                                                                         %
        end                                                                                         %
    %-----------------------------------------------------------------------------------------------%
    case 5                                                                                          %
        while 1                                                                                     %
        choiceHyperPar   = menu('Student-t','alpha_z (default 2.1)', 'beta_z (default 2.1)',...     %
                                            'Haar Level (default 5)','Launch', 'Back');             %
        switch choiceHyperPar                                                                       %
        %-------------------------------------------------------------------------------------------%
        case 1                                                                                      %
            str = '[numerical values must be inputed using the keyboad, then hit enter]';           %
            disp(str)                                                                               %
            azSt0        =  input('Input az0 (default 2.1) :');                                     %
        case 2                                                                                      %
            str = '[numerical values must be inputed using the keyboad, then hit enter]';           %
            disp(str)                                                                               %
            bzSt0        =  input('Input bz0 (default 2.1) :');                                     %
        case 3                                                                                      %
            while 1                                                                                 %
            choiceHaarL  =  menu('Haar Level', '5', '4', '3', 'Back');                              %
            switch choiceHaarL                                                                      %
            case 1                                                                                  %
                l        =  5;                                                                      %
                str=['Haar level is now set at l = ', num2str(l)]; disp(str)                        %
                break;                                                                              %
            case 2                                                                                  %
                l        =  4;                                                                      %
                str=['Haar level is now set at l = ', num2str(l)]; disp(str)                        %
                break;                                                                              %
            case 3                                                                                  %
                l        =  3;                                                                      %
                str=['Haar level is now set at l = ', num2str(l)]; disp(str)                        %
                break;                                                                              %
            case 4                                                                                  %
                break;                                                                              %
            end                                                                                     %
            end                                                                                     %
        case 4                                                                                      %
            %---------------------------------------------------------------------------------------%
            %                                   JMAP Student prior                                  %
            %---------------------------------------------------------------------------------------%
            str=['Bayesian reconstructruction with JMAP estimation and Student prior launched '];   %
            disp(str)                                                                               %
            str=['for data with N = ', num2str(N), ', nproj = ', num2str(nproj),...                 %
                 ' and snr = ', num2str(snr)];                                                      %
            disp(str)                                                                               %
            strPrior     =  '_St_';                                                                 %
            str=['and parameters niter = ', num2str(itermax), ', l = ', num2str(l),...              %
                 ', alpha_z = ', num2str(azSt0), ' and beta_z = ', num2str(bzSt0)];                 %
            disp(str)                                                                               %
            %---------------------------------------------------------------------------------------%
            [MeanTime,fhSt,ve,vxi,vz,V_erf_l2,V_erg_l2,V_erz_l2,V_erf_l1,V_erg_l1,V_erz_l1] ...     %
                         =  TomoGPI_HHBM_St_3d(H_base,g,itermax,fh0_base,f,l,snr,niter_gradient,... %
                                               azSt0,bzSt0);                                        %
            ghSt         =  H_base*fhSt(:); ghSt = reshape(ghSt,H_base.proj_size);                  %
            %---------------------------------------------------------------------------------------%
            figure(300+counter_St); show_obj_3d(fhSt);                                                           %
             RepName=['test_',datestr(now, 'dd-mmm-yyyy-HH-MM')];
            %mkdir([save_path_St,'/',RepName]);%suptitle('Shepp Logan JMAP--St-t reconstructed phantom')  
            %
            save([save_path_St,'/',RepName,'_N', num2str(N), '_nProj', num2str(nproj), '_SNR',...           %
                  num2str(snr), '_nIter', num2str(itermax), '_az', num2str(azSt0), '_bz',...        %
                  num2str(bzSt0), '.mat'],'fhSt','V_erf_l2','V_erg_l2','V_erf_l1','V_erg_l1','ghSt','MeanTime');                                       %
            %---------------------------------------------------------------------------------------%
%             nomf_erf_l2  =  strcat('V_erf_l2_','phantom',int2str(N),strPrior,suffix);               %
%             nomf_erg_l2  =  strcat('V_erg_l2_','phantom',int2str(N),strPrior,suffix);               %
%             nomf_erf_l1  =  strcat('V_erf_l1_','phantom',int2str(N),strPrior,suffix);               %
%             nomf_erg_l1  =  strcat('V_erg_l1_','phantom',int2str(N),strPrior,suffix);               %
%             nomf_fh      =  strcat('fh_','phantom',int2str(N),strPrior,suffix);                     %
%             nomf_gh      =  strcat('gh_','phantom',int2str(N),strPrior,suffix);                     %
%             nomf_time    =  strcat('time_','phantom',int2str(N),strPrior,suffix);                   %
            %---------------------------------------------------------------------------------------%
%             save('-v7.3',[save_path_St,'/',RepName,'/',nomf_erf_l2,'.mat']);                         %
%             save('-v7.3',[save_path_St,'/',RepName,'/',nomf_erg_l2,'.mat']);                         %
%             save('-v7.3',[save_path_St,'/',RepName,'/',nomf_erf_l1,'.mat']);                         %
%             save('-v7.3',[save_path_St,'/',RepName,'/',nomf_erg_l1,'.mat']);                         %
%             save('-v7.3',[save_path_St,'/',RepName,'/',nomf_fh,'.mat'],'fhSt');                                 %
%             save('-v7.3',[save_path_St,'/',RepName,'/',nomf_gh,'.mat']);                                 %
%             save('-v7.3',[save_path_St,'/',RepName,'/',nomf_time,'.mat']);                           %
            %---------------------------------------------------------------------------------------%
            counter_St=counter_St+1;
            break;  
            %
        case 5                                                                                      %
            break;                                                                                  %
        end                                                                                         %
        end                                                                                         %
    %-----------------------------------------------------------------------------------------------%
    case 6                                                                                          %
        while 1                                                                                     %
        choiceHyperPar   = menu('NIG', 'alpha_z (default 2.1)', 'beta_z (default 2.1)',...          %
                                       'Haar Level (default 5)','Launch', 'Back');                  %
        switch choiceHyperPar                                                                       %
        %-------------------------------------------------------------------------------------------%
        case 1                                                                                      %
            str = '[numerical values must be inputed using the keyboad, then hit enter]';           %
            disp(str)                                                                               %
            azNIG0       =  input('Input az0 (default 2.1) :');                                     %
        case 2                                                                                      %
            str = '[numerical values must be inputed using the keyboad, then hit enter]';           %
            disp(str)                                                                               %
            bzNIG0       =  input('Input bz0 (default 2.1) :');                                     %
        case 3                                                                                      %
            while 1                                                                                 %
            choiceHaarL  =  menu('Haar Level', '5', '4', '3', 'Back');                              %
            switch choiceHaarL                                                                      %
            case 1                                                                                  %
                l        =  5;                                                                      %
                str=['Haar level is now set at l = ', num2str(l)]; disp(str)                        %
                break;                                                                              %
            case 2                                                                                  %
                l        =  4;                                                                      %
                str=['Haar level is now set at l = ', num2str(l)]; disp(str)                        %
                break;                                                                              %
            case 3                                                                                  %
                l        =  3;                                                                      %
                str=['Haar level is now set at l = ', num2str(l)]; disp(str)                        %
                break;                                                                              %
            case 4                                                                                  %
                break;                                                                              %
            end                                                                                     %
            end                                                                                     %
        case 4                                                                                      %
            %---------------------------------------------------------------------------------------%
            %                                   JMAP NIG prior                                      %
            %---------------------------------------------------------------------------------------%
            str=['Bayesian reconstructruction with JMAP estimation and NIG prior launched '];       %
            disp(str)                                                                               %
            str=['for data with N = ', num2str(N), ', nproj = ', num2str(nproj),...                 %
                 ' and snr = ', num2str(snr)];                                                      %
            disp(str)                                                                               %
            strPrior     =  '_NIG_';                                                                %
            str=['and parameters niter = ', num2str(itermax), ', l = ', num2str(l),...              %
                 ', alpha_z = ', num2str(azNIG0), ' and beta_z = ', num2str(bzNIG0)];               %
            disp(str)                                                                               %
            %---------------------------------------------------------------------------------------%
            [MeanTime,fhNIG,ve,vxi,vz,V_erf_l2,V_erg_l2,V_erz_l2,V_erf_l1,V_erg_l1,V_erz_l1] ...    %
                         =  TomoGPI_HHBM_NIG_3d(H_base,g,itermax,fh0_base,f,l,snr,niter_gradient,...%
                                                azNIG0,bzNIG0);                                     %
            ghNIG        =  H_base*fhNIG(:); ghNIG  =  reshape(ghNIG,H_base.proj_size);             %
            %---------------------------------------------------------------------------------------%
            figure(400+counter_NIG); show_obj_3d(fhNIG);                                                          %
            suptitle('Shepp Logan JMAP--NIG reconstructed phantom')  
             RepName=['test_',datestr(now, 'dd-mmm-yyyy-HH-MM')];
            %mkdir([save_path_NIG,'/',RepName]);%%
            save([save_path_NIG, '/',RepName, '_NIG_N', num2str(N), '_nProj', num2str(nproj), '_SNR',...         %
                  num2str(snr), '_nIter', num2str(itermax), '_az', num2str(azNIG0), '_bz',...       %
                  num2str(bzNIG0), '.mat'],'fhNIG','V_erf_l2','V_erg_l2','V_erf_l1','V_erg_l1','ghNIG','MeanTime');                                     %
            % save('/espace/tomo_gpi/ResultsMircea/NIG_256_SheppLogan.mat','fhNIG');                %
            %---------------------------------------------------------------------------------------%
%             nomf_erf_l2  =  strcat('V_erf_l2_','phantom',int2str(N),strPrior,suffix);               %
%             nomf_erg_l2  =  strcat('V_erg_l2_','phantom',int2str(N),strPrior,suffix);               %
%             nomf_erf_l1  =  strcat('V_erf_l1_','phantom',int2str(N),strPrior,suffix);               %
%             nomf_erg_l1  =  strcat('V_erg_l1_','phantom',int2str(N),strPrior,suffix);               %
%             nomf_fh      =  strcat('fh_','phantom',int2str(N),strPrior,suffix);                     %
%             nomf_gh      =  strcat('gh_','phantom',int2str(N),strPrior,suffix);                     %
%             nomf_time    =  strcat('time_','phantom',int2str(N),strPrior,suffix);                   %
%             %---------------------------------------------------------------------------------------%
%             save('-v7.3',[save_path_NIG,'/',RepName,'/',nomf_erf_l2,'.mat'],'V_erf_l2');                        %
%             save('-v7.3',[save_path_NIG,'/',RepName,'/',nomf_erg_l2,'.mat'],'V_erg_l2');                        %
%             save('-v7.3',[save_path_NIG,'/',RepName,'/',nomf_erf_l1,'.mat'],'V_erf_l1');                        %
%             save('-v7.3',[save_path_NIG,'/',RepName,'/',nomf_erg_l1,'.mat'],'V_erg_l1');                        %
%             save('-v7.3',[save_path_NIG,'/',RepName,'/',nomf_fh,'.mat'],'fhNIG');                               %
%             save('-v7.3',[save_path_NIG,'/',RepName,'/',nomf_gh,'.mat'],'ghNIG');                               %
%             save('-v7.3',[save_path_NIG,'/',RepName,'/',nomf_time,'.mat'],'MeanTime');                          %
            %---------------------------------------------------------------------------------------%
           counter_NIG=counter_NIG+1;
            break;                                                                                  %
        case 5                                                                                      %
            break;                                                                                  %
        end                                                                                         %
        end                                                                                         %
    %-----------------------------------------------------------------------------------------------%
    case 7                                                                                          %
        while 1                                                                                     %
        choiceHyperPar   = menu('VG', 'alpha_z (default 2.1)', 'beta_z (default 1/2.1)',...         %
                                      'Haar Level (default 5)','Launch','Break');                   %
        switch choiceHyperPar                                                                       %
        %-------------------------------------------------------------------------------------------%
        case 1                                                                                      %
            str = '[numerical values must be inputed using the keyboad, then hit enter]';           %
            disp(str)                                                                               %
            azVG0        = input('Input az0 (default 2.1) :');                                      %
        case 2                                                                                      %
            str = '[numerical values must be inputed using the keyboad, then hit enter]';           %
            disp(str)                                                                               %
            bzVG0        = input('Input bz0 (default 1/2.1) :');                                    %
        case 3                                                                                      %
            while 1                                                                                 %
            choiceHaarL  =  menu('Haar Level', '5', '4', '3', 'Back');                              %
            switch choiceHaarL                                                                      %
            case 1                                                                                  %
                l        =  5;                                                                      %
                str=['Haar level is now set at l = ', num2str(l)]; disp(str)                        %
                break;                                                                              %
            case 2                                                                                  %
                l        =  4;                                                                      %
                str=['Haar level is now set at l = ', num2str(l)]; disp(str)                        %
                break;                                                                              %
            case 3                                                                                  %
                l        =  3;                                                                      %
                str=['Haar level is now set at l = ', num2str(l)]; disp(str)                        %
                break;                                                                              %
            case 4                                                                                  %
                break;                                                                              %
            end                                                                                     %
            end                                                                                     %
        case 4                                                                                      %
            %---------------------------------------------------------------------------------------%
            %                                   JMAP VG prior                                       %
            %---------------------------------------------------------------------------------------%
            str=['Bayesian reconstructruction with JMAP estimation and VG prior launched '];        %
            disp(str)                                                                               %
            str=['for data with N = ', num2str(N), ', nproj = ', num2str(nproj),...                 %
                 ' and snr = ', num2str(snr)];                                                      %
            disp(str)                                                                               %
            strPrior     =  '_VG_';                                                                 %
            str=['and parameters niter = ', num2str(itermax), ', l = ', num2str(l),...              %
                 ', alpha_z = ', num2str(azVG0), 'and beta_z = ', num2str(bzVG0)];                  %
            disp(str)                                                                               %
            %---------------------------------------------------------------------------------------%
            [MeanTime,fhVG,ve,vxi,vz,V_erf_l2,V_erg_l2,V_erz_l2,V_erf_l1,V_erg_l1,V_erz_l1] ...     %
                         =  TomoGPI_HHBM_VG_3d(H_base,g,itermax,fh0_base,f,l,snr,niter_gradient,... %
                                               azVG0,bzVG0);                                        %
            ghVG         =  H_base*fhVG(:); ghVG  =  reshape(ghVG,H_base.proj_size);                %
            %---------------------------------------------------------------------------------------%
            figure(500+counter_VG); show_obj_3d(fhVG);                                                           %
            suptitle('Shepp Logan JMAP--VG reconstructed phantom')    
             RepName=['test_',datestr(now, 'dd-mmm-yyyy-HH-MM')];
            %mkdir([save_path_VG,'/',RepName]);%%
            save([save_path_VG,'/',RepName, '_VG_N', num2str(N), '_nProj', num2str(nproj), '_SNR',...           %
                  num2str(snr), '_nIter', num2str(itermax), '_az', num2str(azNIG0), '_bz',...       %
                  num2str(bzVG0), '.mat'],'fhVG','V_erf_l2','V_erg_l2','V_erf_l1','V_erg_l1','ghVG','MeanTime');                                       %
            % save('/espace/tomo_gpi/ResultsMircea/VG_256_SheppLogan.mat','fhVG');                  %
            %---------------------------------------------------------------------------------------%
%             nomf_erf_l2  =  strcat('V_erf_l2_','phantom',int2str(N),strPrior,suffix);               %
%             nomf_erg_l2  =  strcat('V_erg_l2_','phantom',int2str(N),strPrior,suffix);               %
%             nomf_erf_l1  =  strcat('V_erf_l1_','phantom',int2str(N),strPrior,suffix);               %
%             nomf_erg_l1  =  strcat('V_erg_l1_','phantom',int2str(N),strPrior,suffix);               %
%             nomf_fh      =  strcat('fh_','phantom',int2str(N),strPrior,suffix);                     %
%             nomf_gh      =  strcat('gh_','phantom',int2str(N),strPrior,suffix);                     %
%             nomf_time    =  strcat('time_','phantom',int2str(N),strPrior,suffix);                   %
%             %---------------------------------------------------------------------------------------%
%             save('-v7.3',[save_path_VG,'/',RepName,'/',nomf_erf_l2,'.mat'],'V_erf_l2');                         %
%             save('-v7.3',[save_path_VG,'/',RepName,'/',nomf_erg_l2,'.mat'],'V_erg_l2');                         %
%             save('-v7.3',[save_path_VG,'/',RepName,'/',nomf_erf_l1,'.mat'],'V_erf_l1');                         %
%             save('-v7.3',[save_path_VG,'/',RepName,'/',nomf_erg_l1,'.mat'],'V_erg_l1');                         %
%             save('-v7.3',[save_path_VG,'/',RepName,'/',nomf_fh,'.mat'],'fhVG');                                 %
%             save('-v7.3',[save_path_VG,'/',RepName,'/',nomf_gh,'.mat'],'ghVG');                                 %
%             save('-v7.3',[save_path_VG,'/',RepName,'/',nomf_time,'.mat'],'MeanTime');                           %
            %---------------------------------------------------------------------------------------%
              counter_VG=counter_VG+1;
            break;                                                                                  %
        case 5                                                                                      %
            break;                                                                                  %
        end                                                                                         %
        end                                                                                         %
    case 8                                                                                          %
        break;                                                                                    %
    end                                                                                             %
    end                                                                                             %
%---------------------------------------------------------------------------------------------------%
%                                                                                                   %
case 4                                                                                              %
    close all;clear all;break;
    
end                                                                                                 %
end                                                                                                 %
%---------------------------------------------------------------------------------------------------%