clc; clear; close all

addpath("D:\My_Study\HuaLien_April3_M7_4\Prog-1\lib\gmm");
    
addpath("D:\My_Study\HuaLien_April3_M7_4\References");

df = readtable("April12M7_4hGM_CWBTSMIP.csv");

M0 =  7.0;

ind = logical((isnan(df.VS30)) +  (df.VS30 < 0));

df.VS30(ind) = 500;

df.Z1_0(ind) = exp(-3.73/2*log((df.VS30(ind).^2+290.53^2)/(1750^2+290.53^2)));

inp = df.Z1_0 < 0;

df.Z1_0(inp) = exp(-3.73/2*log((df.VS30(inp).^2+290.53^2)/(1750^2+290.53^2)));

Ztor = (max(3.5384-2.60*max(M0-5.8530,0),0)).^2;
%% STAID of pulse like Directivity and Pulse like motion 


%% Sa-tables
colnames = df.Properties.VariableNames;
colnames_period = colnames(startsWith(colnames,'T'));


VarNames = df(:,colnames_period).Properties.VariableNames;
period_loc = regexp(VarNames,'\d*','Match');
% Obtain the frequencies range.
period = nan(1, length(period_loc));
 
for k = 1:length(period_loc)
    if numel(period_loc{1,k})>1
       period(1,k) = str2double(strcat(period_loc{1,k}(1),'.',period_loc{1,k}(2)));
    else
       period(1,k) = str2double(period_loc{1,k});
    end
end
%% 
sd1 = [121.37 23.87 2.85];
sd2 = [121.61 24.32 2.85];
sd3 = [121.87 24.21 48.95];
sd4 = [121.62 23.76 48.95];
sd5 = [121.37 23.87 2.85];

Flat = [24.382, 24.221, 23.596, 23.758, 24.382];
Flon = [121.638, 121.999, 121.669, 121.306, 121.638];


% figure;
% geoscatter(df.StaLat,df.StaLon, 10, df.T0_01s,'filled'); hold on;
% geoplot(Flat, Flon,'r-'); hold on;
% geoplot(df.EveLat(1), df.EveLon(1),'pr','MarkerSize',12);


%%
figure;
for M0 = [7.4, 7.0]
    dfc = df(df.MW == M0, :);
    Satable = table2array(dfc(:,colnames_period));

    dip = 45;
    lambda = 90;
    sigma_obs = [];
    sigma_pred = [];
    
    for T = period
        yprd = nan(numel(dfc.MW), 1);
        sigma_ask14 = nan(numel(dfc.MW), 1);
        
        for i = 1:length(dfc.MW)
            [yprd(i,1), sigma_ask14(i,1)] = ASK_2014_nga(dfc.MW(i), T, dfc.Rrup(i), 0, 0, 0, dfc.Ztor(i), dip, lambda, 0, 0, 999, dfc.Z1_0(i), dfc.VS30(i), 0, 0);
        end
    
        Resid = log(Satable(:, period==T)./980.6) - log(yprd);
        sigma_obs = cat(1, sigma_obs, std(Resid));
        sigma_pred = cat(1, sigma_pred, mean(sigma_ask14));
    end

    if M0 == 7.4
        a=semilogx(period,sigma_obs,'b','LineWidth',2); hold on;
        b=semilogx(period,sigma_pred,'r','LineWidth',2); 
    else
        c=semilogx(period,sigma_obs,'b--','LineWidth',2); hold on;
        d=semilogx(period,sigma_pred,'r--','LineWidth',2); 
    end

end

legend([a,b,c,d], 'data, M = 7.4','predicted, M = 7.4',...
                  'data, M = 7.0','predicted, M = 7.0');
xlabel('Period (s)','FontSize',12); grid on;
ylabel('Standard deviation, \sigma','FontSize',15);


%% Sigma Magnitude - dependent of ASK14. 

fig = figure;
for M = [4, 5, 6, 7, 8]
        sigmaM = [];
    for T = period
        [~, sigma_ask] = ASK_2014_nga(M, T, 30, 0, 0, 0, 999, dip, lambda, 0, 0, 999, 999, 760, 0, 0);
        sigmaM = cat(1,sigmaM,sigma_ask);
    end
    k = find([4, 5, 6, 7, 8] ==M);
    a(k) = semilogx(period,sigmaM,'LineWidth',2); hold on;
end
legend(a, 'M = 4.0','M = 5.0',...
                  'M = 6.0', 'M = 7.0','M = 8.0');
xlabel('Period (s)','FontSize',12); grid on;
ylabel('Standard deviation, \sigma','FontSize',15);

axis([0.01, 10, 0.3, 1.2]);

fig_name = 'Sigma_Model';
saveas(fig,[pwd,'\plots_h2v_ila\',fig_name,'.png']);