close all; clear all; clc;

%%% Generate the PGA values for Kontum earthquake using the typical GPMEs
%%% for Vietnam (following Phung et al., 2024)
addpath('.\src');
%%% make output directory
dirPath="./gmpes"
if ~exist(dirPath, 'dir')
    % If the directory does not exist, create it
    mkdir(dirPath);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ------------- Pre-condition several parameters ------------------------
%% https://earthquake.usgs.gov/earthquakes/eventpage/us6000ngfw/moment-tensor
Vs30=750; %km/s hard rock
dist=0:1:300; % km
Ztor=10;
delta=58; % average dip of the rupture place (degree)
lambda=-79; % rake angle (degree)
arb=0; % average
Rrup=10;
Rjb=10;
Zvs=999; % No sediment effect
A1100=0;
period=0; %PGA
% ----------- ask 14
region=0; %global
fas= 0; % Flag for aftershocks
HW = 0; % Falg for hanging wall sites
Z10=Zvs; % Basin depth (km);
W=1;% Down-dip rupture width (km)
FVS30 = 0;
% --------------- %% ---------------
option=3;
Nstd=0;
%%% 
% % magnitude range generation
M=5.0:0.01:8.5;
% buiding empty array
for M_now = M
    Sa_M=[];
    Sigma_M=[];
    for dist_now = dist
        
        Rrup=dist_now;
        Rjb=dist_now;
        % [Sa_cb, sigma, period1] = CB_2008_nga (M_now, period, Rrup, Rjb, Ztor, delta, lambda, Vs30, Zvs, arb, A1100);
        
         [Sa_cb, sigma] = ASB_2014(M_now,period,Rjb,Vs30,lambda,option,Nstd); 
         
%         [Sa_ask] = ASK_2014_nga(M, period, Rrup, Rjb, 0, Rrup, Ztor, delta, lambda, fas, HW, W, Z10, Vs30, FVS30, region)
        % convert from G to gal
        Sa_cb = Sa_cb*98;
%         Sa_ask = Sa_ask*98;
        sigma = sigma;
        % append the value into the array
        Sa_M = [Sa_M,Sa_cb]; % Sa in gal
        Sigma_M = [Sigma_M,sigma];
    end 

    % Combine the arrays into a matrix with two columns
    data_CB08 = [dist',Sa_M', Sigma_M'];
    filename_CB08=strjoin([dirPath,"/CB08_M", num2str(M_now, '%.2f'), "_PGA.dat"],"")
%     filename_ASK14=strjoin([dirPath,"/ASK14_M", num2str(M_now, '%.2f'), "_PGA.dat"],"")

    disp(filename_CB08)
    % Write the combined data to an ASCII file using dlmwrite
    dlmwrite(filename_CB08, data_CB08, 'delimiter', ' ');
    clear data_CB08
end
% [Sa, sigma, period1] = CB_2008_nga (M_now, 0, Rrup, Rjb, Ztor, delta, lambda, 750, Zvs, arb,A1100)