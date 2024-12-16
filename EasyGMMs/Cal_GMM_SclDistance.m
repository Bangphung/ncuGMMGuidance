clc; clear; close all; 

addpath("D:\My_Study\HuaLien_April3_M7_4\Prog-1\lib\gmm");

%% Save data
    file_dir = fileparts(mfilename ('fullpath'));
    if ~isempty(file_dir)
        file_dir = [file_dir,'\'];
    end
    dir_output = [file_dir,sprintf('plots_gmm/')];
    
    if ~(exist(dir_output, 'dir')==7)
        mkdir(dir_output)
    end 

%% Given Scenarios as input for the GMMs
M = 7;
Ztor = 0; 
Vs30 = 760;
dist_range = logspace(0,3,100);

delta = 45;
lambda = 90;
VS30 = 760;
T = 0;
region = 0;
FVS30 = 0;

Rjb = 999; 
Rx = 999;
Ry0 = 999;
fas = 0;
HW = 0;
W = 999;
Z1_0 = 999;
Zbot = 999;
Z25 = 999;
Zhyp = 999;

% Input for Ch20
FRO = 1; 
FSS = 0;
FNO = 0; 
Finter = 0; 
Fintra = 0; 
Fas = 0; 
Fma = 0;  
 color = lines(8);
%% Create vector to store the predicted values
    Sa_ask = nan(length(dist_range),1);
    Sa_bssa = nan(length(dist_range),1);
    
    Sa_cb = nan(length(dist_range),1);
    Sa_cy = nan(length(dist_range),1);
    
    Sa_ch = nan(length(dist_range),1);
    Sa_ph = nan(length(dist_range),1);
    Sa_cb08 = nan(length(dist_range),1);
    Sa_asb14 = nan(length(dist_range),1);
%% Periods considered for the predictions
periods = [0.01 0.2 1];

fig = figure('position',[50 50 1400 400]);
t = tiledlayout(1,3,"TileSpacing","compact");
for k = 1:3
    nexttile
    T = periods(k);
    
    for i = 1:length(dist_range)
        Sa_ask(i,1) = ASK_2014_nga(M, T, dist_range(i), Rjb, Rx, Ry0, 999, delta, lambda, fas, HW, W, Z1_0, VS30, FVS30, region);
        Sa_bssa(i,1) = BSSA_2014_nga(M, T, dist_range(i), lambda, region, Z1_0, VS30);
        Sa_cb(i,1) = CB_2014_nga(M, T, dist_range(i), Rjb, Rx, W, Ztor, Zbot, delta, lambda, HW, VS30, Z25, Zhyp, region);
        Sa_cy(i,1) =  CY_2014_nga(M, T, dist_range(i), Rjb, Rx, Ztor, delta, lambda, Z1_0, VS30, HW, FVS30, region);
        Sa_ch(i,1) = Chao19H(T,M,dist_range(i), Ztor, VS30,Z1_0, FRO,FSS,FNO,Finter,Fintra,Fas,Fma,FVS30);
        Sa_ph(i,1) = Phung_2018c2_NGAw2_TW(M, T, dist_range(i), 0, 0, Ztor, 90, 0, 0, Vs30, Z1_0, 1, region);
        
        option = 1; % for Rjb model;
        Nstd = 1; % for number of std;
        Sa_asb14(i,1) = ASB_2014(M, T, dist_range(i), Vs30, lambda,option,Nstd);
    end
    
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    loglog(dist_range,Sa_ask,'b-','LineWidth',2); hold on;
    loglog(dist_range,Sa_bssa,'c-','LineWidth',2); hold on;
    loglog(dist_range,Sa_cb,'g-','LineWidth',2); hold on;
    loglog(dist_range,Sa_cy,'m-','LineWidth',2); hold on;
    loglog(dist_range,Sa_ch,'k-','LineWidth',2); hold on;
    loglog(dist_range,Sa_ph,'r','LineWidth',2); hold on;
    loglog(dist_range,Sa_asb14,'b--','LineWidth',2); hold on;

    xlabel('R_{rup}, [km]');
    if k==1
        ylabel('PGA, [gal]')
    else
        ylabel(['5% damped PSA at ',num2str(T),'s, [gal]']);
    end
    
    title('SS: M=7.0; Z_{tor}=0 km; V_{S30}=760m/s');
    axis([1 1000 1e-4 3]);
    grid on;
    
    if k==1
        legend('ASK14','BSSA14','CB14','CY14','CH20','PH20','ASB14',...
               'location','best','NumColumns', 3);
    else
        legend off
    end

end

fig_name = 'DistanceScalingGMMs.png';
saveas(fig, [dir_output, fig_name]);


%% Magnitude Scaling 
mag_range = 4:0.1:8;
Rrup = 30; 
%% Create space to store the predicitons
   Sa_ask = nan(length(mag_range),1);
   Sa_bssa = nan(length(mag_range),1);
    
   Sa_cb = nan(length(mag_range),1);
   Sa_cy = nan(length(mag_range),1);
    
   Sa_ch = nan(length(mag_range),1);
   Sa_ph = nan(length(mag_range),1);
   Sa_cb08 = nan(length(mag_range),1);
   Sa_asb14 = nan(length(mag_range),1);
%% Plot the Magnitude Scaling 
fig = figure('position',[50 50 1400 400]);
t = tiledlayout(1,3,"TileSpacing","compact");
for k = 1:3
    nexttile
    T = periods(k);
    for i = 1:length(mag_range)
        % Ztor computed from M-Ztor relationship of CY14
        Ztor = (max(2.704-1.226*max(mag_range(i)-5.849,0),0))^2;
        Rjb = sqrt(Rrup.^2 - Ztor.^2);
        Sa_ask(i,1) = ASK_2014_nga(mag_range(i), T, Rrup, Rjb, Rx, Ry0, 999, delta, lambda, fas, HW, W, Z1_0, VS30, FVS30, region);
        Sa_bssa(i,1) = BSSA_2014_nga(mag_range(i), T, Rrup, lambda, region, Z1_0, VS30);
        Sa_cb(i,1) = CB_2014_nga(mag_range(i), T, Rrup, Rjb, Rx, W, Ztor, Zbot, delta, lambda, HW, VS30, Z25, Zhyp, region);
        Sa_cy(i,1) =  CY_2014_nga(mag_range(i), T, Rrup, Rjb, Rx, Ztor, delta, lambda, Z1_0, VS30, HW, FVS30, region);
        Sa_ch(i,1) = Chao19H(T, mag_range(i), Rrup, Ztor, VS30,Z1_0, FRO,FSS,FNO,Finter,Fintra,Fas,Fma,FVS30);
        Sa_ph(i,1) = Phung_2018c2_NGAw2_TW(mag_range(i), T, Rrup, 0, 0, Ztor, 90, 0, 0, Vs30, Z1_0, 1, region);
        
        option = 1; % for Rjb model;
        Nstd = 1; % for number of std;
        Sa_asb14(i,1) = ASB_2014(mag_range(i), T, Rjb, Vs30, lambda,option,Nstd);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    semilogy(mag_range,Sa_ask,'b-','LineWidth',2); hold on;
    semilogy(mag_range,Sa_bssa,'c-','LineWidth',2); hold on;
    semilogy(mag_range,Sa_cb,'g-','LineWidth',2); hold on;
    semilogy(mag_range,Sa_cy,'m-','LineWidth',2); hold on;
    semilogy(mag_range,Sa_ch,'k-','LineWidth',2); hold on;
    semilogy(mag_range,Sa_ph,'r','LineWidth',2); hold on;
    semilogy(mag_range,Sa_asb14,'b--','LineWidth',2); hold on;

    xlabel('Magnitude, [M]');
    if k==1
        ylabel('PGA, [gal]')
    else
        ylabel(['5% damped PSA at ',num2str(T),'s, [gal]']);
    end
    
    title('SS: R_{rup}=30; Z_{tor} = f_{CY14}(M) km; V_{S30}=760m/s');
    axis([4 8 3*1e-3 3]);
    grid on;
    
    if k==1
        legend('ASK14','BSSA14','CB14','CY14','CH20','PH20','ASB14',...
               'location','best','NumColumns', 3);
    else
        legend off
    end

end

fig_name = 'MagScalingGMMs.png';
saveas(fig, [dir_output, fig_name]);
