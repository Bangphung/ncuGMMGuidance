clc; clear; close all;

addpath("D:\My_Study\HuaLien_April3_M7_4\Prog-1\lib\gmm");

addpath("D:\My_Study\HuaLien_April3_M7_4\AnaTaipeiSiteResponse");

df1 = readtable("20180206HL_hGM_CWBTSMIP.csv");

df2 = readtable("ChihShangM7_hGM_CWBTSMIP.csv");

df3 = readtable("April12M7_4hGM_CWBTSMIP.csv");

df3 = df3(df3.MW~=7, :);

df1E = df1(:,{'EQID','EveLat','EveLon','STAID','StaLat','StaLon','MW','Depth','Rhyp','VS30','Z1_0'});
df1E = renamevars(df1E,"Rhyp","Rrup");
df2E = df2(:,{'EQID','EveLat','EveLon','STAID','StaLat','StaLon','MW','Depth','Rhyp','VS30','Z1_0'});
df2E = renamevars(df2E,"Rhyp","Rrup");

df3E = df3(:,{'EQID','EveLat','EveLon','STAID','StaLat','StaLon','MW','Depth','Rrup','VS30','Z1_0'});

df1psa = table2array(df1(:,19:end))./980.6;

df2psa = table2array(df2(:,22:end))./980.6;

df3psa = table2array(df3(:,22:end))./980.6;
%% Save data
    file_dir = fileparts(mfilename ('fullpath'));
    if ~isempty(file_dir)
        file_dir = [file_dir,'\'];
    end
    dir_output = [file_dir,sprintf('plots_gmm/')];
    
    if ~(exist(dir_output, 'dir')==7)
        mkdir(dir_output)
    end 
%% TSSHAC level3 dataset.
load('SSHAC_GM_v9.mat');

ds = dataset2table(ds);

id = ds.MW >= 5.5 & ds.Rrup <= str2double(ds.Rmax1);  

ip = startsWith(ds.eq_type,'shallow crustal'); % shallow

ipd = logical(id.*ip);

dsE = ds(ipd,{'EQSN', 'Hyp_Lat','Hyp_Long','STA_ID','STA_Lat_Y', 'STA_Lon_X', 'MW', 'Ztor','Rrup','Vs30','Z1_0'});

dsE = renamevars(dsE,["EQSN", "Hyp_Lat","Hyp_Long","STA_ID","STA_Lat_Y", "STA_Lon_X", "Ztor", "Vs30"],...
                      ["EQID", "EveLat", "EveLon", "STAID","StaLat","StaLon","Depth","VS30"]);

dsE.Z1_0 = str2double(dsE.Z1_0);

dSpsa = table2array(ds(ipd, 80:184));

dpsa = [df1psa; df2psa; df3psa; dSpsa];

dpsa = dpsa(dpsa(:,1)>0,:);

deq = [df1E; df2E; df3E; dsE];

deq = deq(dpsa(:,1)>0,:);

deq.VS30(isnan(deq.VS30)) = 500;

deq.VS30(deq.VS30 <0 ) = 500;
%% Plot the Events - Map 
% [EQID, eqid] = unique(deq.EQID);
% 
% [STAID, staid] = unique(deq.STAID);
% 
% 
% fig1 = figure;
% a = geoplot(deq.StaLat(staid), deq.StaLon(staid),'k^','MarkerFaceColor',[0.5, 0.5, 0.5],'MarkerSize',1); hold on;
% b = geoscatter(deq.EveLat(eqid), deq.EveLon(eqid), 10*deq.MW(eqid), deq.MW(eqid),'o','filled'); hold on;
% legend([a, b],'stations','Events (M_{W}=5.5-7.6)','Location','best');
% c = colorbar;
% c.Label.String = "Magnitude";
% 
% geobasemap topographic
% 
% fig_name = 'LArgeEQDistibution.png';
% saveas(fig1, [dir_output, fig_name]);
%% Plot M-R distribution
% fig2 = figure;
% semilogx(deq.Rrup, deq.MW, 'ko','MarkerFaceColor', 'b');
% xlabel('Distance (km)');
% ylabel('Magnitude, M_W');
% axis([0.1 500, 5.4, 8]);
% title('Taiwan EQ, M_W = 5.5 - 7.6');
% grid on;
% 
% fig_name = 'MagDistanceDistibution.png';
% saveas(fig2, [dir_output, fig_name]);
% 
% fig3 = figure;
% plot(deq.Depth, deq.MW, 'ko','MarkerFaceColor', 'b','MarkerSize',10);
% xlabel('Depth, Z_{tor} (km)');
% ylabel('Magnitude, M_W');
% axis([-0.5, 35, 5.4, 8]);
% title('Taiwan EQ, M_W = 5.5 - 7.6');
% grid on;
% fig_name = 'MagDepthDistibution.png';
% saveas(fig3, [dir_output, fig_name]);


%% Perform Residual Analysis

HW = 0; 
Rjb = 0;
Rx = 0;
Ry0 = 0;
W = 0;
delta = 90;
lambda = 0;
fas = 0;
FVS30 = 0;
region = 1; % for Taiwan
T = 0; % for PGA

% Given observed data 
  ydata = dpsa(:,1); % PGA

  model = 'Ph20';

if strcmp(model,'Ph20')
    Sa = nan(length(deq.EQID),1);
    for i = 1:length(deq.EQID)
        Sa(i,1) = Phung_2019h_NGAw2_TW(deq.MW(i), T, deq.Rrup(i), Rjb, Rx, deq.Depth(i), delta, HW, lambda, deq.VS30(i),-999, region, fas);
    end

else
    Sa = nan(length(deq.EQID),1);
    for i = 1:length(deq.EQID)
        Sa(i,1) = ASK_2014_nga(deq.MW(i), T, deq.Rrup(i), Rjb, Rx, Ry0, deq.Depth(i), delta, lambda, fas, HW, W, 999, deq.VS30(i), FVS30, 0);
    end
end
% Calculate the total Residuals.


Resid = log(ydata) - log(Sa);

% Create table to perform the mixeffs Regression.

tab = table(Resid, deq.EQID, deq.STAID,'VariableNames',{'Resid','EQID','STAID'});

lme = fitlme(tab, 'Resid ~ 1 + (1|EQID) + (1|STAID)');

 [randeff, Bnames] = randomEffects(lme); 

 dBe = randeff(ismember(Bnames.Group,'EQID'));

 dS2S = randeff(ismember(Bnames.Group,'STAID'));

 dWes = residuals(lme); 

 c0 = fixedEffects(lme);

 [EQID, ia] = unique(deq.EQID); Mag = deq.MW(ia); Zh = deq.Depth(ia);
              
 [STAID, ib] = unique(deq.STAID); VS30 = deq.VS30(ib); Z1_0 = deq.Z1_0(ib);            


%% Plot the residuals versus M, Rrup and VS30
    meanfunc = @(x,y)[mean(x), mean(y), 1.96*std(y)./numel(y)];
    
    igpm = discretize(Mag, [5.45, 6.0, 6.5, 7.0, 7.7]);
    
    ermag = splitapply(meanfunc,Mag,dBe,igpm);
    
    igpZ = discretize(Zh,[0, 5, 10, 15, 20, 25, 35]);
    
    erZh = splitapply(meanfunc,Zh,dBe,igpZ);
    
    igpR = discretize(deq.Rrup, cat(2,logspace(-1, 2.0, 7), 150, 250));

    erR = splitapply(meanfunc,deq.Rrup,dWes,igpR);

%% Perform Linear Fit to the residuals 
     Xnew = (5.4:0.1:8)';
    [ypred, yci] = fitLN(Mag,dBe,Xnew);
    
     X2new = (0:0.5:35)';
    [y2pred, y2ci] = fitLN(Zh,dBe,X2new);
    
     X3new = (logspace(-2,3,100))';
    [y3pred, y3ci] = fitLN(deq.Rrup, dWes, X3new);

%% Plot the residuals
figure('Position',[20, 20, 600, 400]);
plot(Mag, dBe, 'ko','MarkerFaceColor','b'); hold on;
e = errorbar(ermag(:,1),ermag(:,2),ermag(:,3),'sqr');
e.CapSize = 12; e.LineWidth = 1.2; e.MarkerSize = 10;
plot(Xnew,ypred,'r','LineWidth',1.5); hold on;
plot(Xnew,yci,'r--','LineWidth',1.5); 
xlabel('Magnitude');
ylabel('PGA, \delta Be');
axis([5.4, 8, -2.0, 2.0]);
title(model);
grid on;



figure('Position',[20, 20, 600, 400]);
plot(Zh, dBe, 'ko','MarkerFaceColor','b'); hold on;
e = errorbar(erZh(:,1),erZh(:,2),erZh(:,3),'sqr');
e.CapSize=12; e.LineWidth = 1.2; e.MarkerSize = 10;
plot(X2new,y2pred,'r','LineWidth',1.5); hold on;
plot(X2new,y2ci,'r--','LineWidth',1.5); 
xlabel('Depth, Z_{TOR} (km)');
ylabel('PGA, \delta Be');
axis([0, 35, -2.0, 2.0]);
title(model);
grid on;


figure('Position',[20, 20, 600, 400]);
semilogx(deq.Rrup, dWes, 'ko','MarkerFaceColor','b'); hold on;
e = errorbar(erR(:,1),erR(:,2),erR(:,3),'sqr');
e.CapSize = 12; e.LineWidth = 1.2; e.MarkerSize = 10;
plot(X3new,y3pred,'r','LineWidth',1.5); hold on;
plot(X3new,y3ci,'r--','LineWidth',1.5); 
xlabel('Rrupture distance R_{rup}, (km)');
ylabel('PGA, \delta WS_{es}');
axis([1, 500, -4.5, 4.5]);
title(model);
grid on;



%% Save data. 
ResMat = ones(length(deq.EQID),3);
for j = 1:length(EQID)
    ind = ismember(deq.EQID, EQID(j));
    ResMat(ind,1) = dBe(j);
end

for k = 1:length(STAID)
    ink = ismember(deq.STAID, STAID(k));
    ResMat(ink,2) = dS2S(k);
end

ResMat(:,3) = dWes;

ResMat2 = array2table(ResMat,'VariableNames',{'dBe','dS2S','dWSes'});

deq.SSN = findgroups(deq.STAID);

dres = [deq(:,{'EQID','SSN','MW','Depth','Rrup','VS30'}) ResMat2];

% writetable(dres,'PGAResMatPh20.csv');


%% Create funstion
function [ypred,yci1] = fitLN(X,y,Xnew)
     mdl = fitlm(X, y); 
     ypred = predict(mdl,Xnew);
     [~,yci1] = predict(mdl,Xnew,'Alpha',0.1,'Simultaneous',true);
end