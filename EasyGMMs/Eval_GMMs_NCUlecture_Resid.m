clc; clear; close all;

    model = 'Ph20, PGA';

if strcmp(model,'Ph20, PGA')
    df = readtable("PGAResMatPh20.csv");
else
    df = readtable("PGAResMatASK14.csv");
end

     [EQID, ia] = unique(df.EQID); Mag = df.MW(ia); Zh = df.Depth(ia);
                  
     [STAID, ib] = unique(df.SSN); VS30 = df.VS30(ib); Rrup = df.Rrup;

     dBe = df.dBe(ia); dS2S = df.dS2S(ib); dWSes = df.dWSes;
%% Plot the residuals versus M, Rrup and VS30
    meanfunc = @(x,y)[mean(x), mean(y), 1.96*std(y)./numel(y)];
    
    igpm = discretize(Mag, [5.45, 6.0, 6.5, 7.0, 7.7]);
    
    ermag = splitapply(meanfunc,Mag,dBe,igpm);
    
    igpZ = discretize(Zh,[0, 5, 10, 15, 20, 25, 35]);
    
    erZh = splitapply(meanfunc,Zh,dBe,igpZ);
    
    igpR = discretize(Rrup, cat(2,logspace(-1, 2.0, 7), 150, 250));

    erR = splitapply(meanfunc,Rrup,df.dWSes,igpR);


    igpV = discretize(VS30, cat(2, logspace(2, 3.15, 8),  1550));

    erV = splitapply(meanfunc,VS30,dS2S,igpV);
%% Perform Linear Fit to the residuals 
     Xnew = (5.4:0.1:8)';
    [ypred, yci] = fitLN(Mag,dBe,Xnew);
    
     X2new = (0:0.5:35)';
    [y2pred, y2ci] = fitLN(Zh,dBe,X2new);
    
     X3new = (logspace(-2,3,100))';
    [y3pred, y3ci] = fitLN(Rrup, dWSes, X3new);

     X4new = (logspace(2,3.25,100))';
    [y4pred, y4ci] = fitLN(VS30, dS2S, X4new);

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
semilogx(Rrup, dWSes, 'ko','MarkerFaceColor','b'); hold on;
e = errorbar(erR(:,1),erR(:,2),erR(:,3),'sqr');
e.CapSize = 12; e.LineWidth = 1.2; e.MarkerSize = 10;
plot(X3new,y3pred,'r','LineWidth',1.5); hold on;
plot(X3new,y3ci,'r--','LineWidth',1.5); 
xlabel('Rrupture distance R_{rup}, (km)');
ylabel('PGA, \delta WS_{es}');
axis([1, 500, -4.5, 4.5]);
title(model);
grid on;



figure('Position',[20, 20, 600, 400]);
semilogx(VS30, dS2S, 'ko','MarkerFaceColor','b'); hold on;
e = errorbar(erV(:,1), erV(:,2), erV(:,3), 'sqr');
e.CapSize = 12; e.LineWidth = 1.2; e.MarkerSize = 10;
plot(X4new,y4pred,'r','LineWidth',1.5); hold on;
plot(X4new,y4ci,'r--','LineWidth',1.5); 
xlabel('V_{S30}, (m/s)');
ylabel('PGA, \delta S2S_{s}');
axis([100, 1500, -2.5, 2.5]);
title(model);
grid on;

%% Create funstion
function [ypred,yci1] = fitLN(X,y,Xnew)
     mdl = fitlm(X, y); 
     ypred = predict(mdl,Xnew);
     [~,yci1] = predict(mdl,Xnew,'Alpha',0.1,'Simultaneous',true);
end