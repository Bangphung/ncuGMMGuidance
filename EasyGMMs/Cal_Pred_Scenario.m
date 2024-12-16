clc; clear; close all
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
%%    
M = 7.0; 
Rrup = 30; 
Ztor = 0;
delta = 45;
lambda = 90;
VS30 = 760;
T = 0;
region = 0;
FVS30 = 0;

Rjb = sqrt(Rrup.^2 - Ztor.^2); 
Rx = 999;
Ry0 = 999;
fas = 0;
HW = 0;
W = 999;
Z1_0 = 999;
Zbot = 999;
Z25 = 999;
Zhyp = 999;

[SaASK14, sig_ask14] = ASK_2014_nga(M, T, Rrup, Rjb, Rx, Ry0, 999, delta, lambda, fas, HW, W, Z1_0, VS30, FVS30, region);
[SaCY14, sig_cy14] = CY_2014_nga(M, T, Rrup, Rjb, Rx, Ztor, delta, lambda, Z1_0, VS30, HW, FVS30, region);
[SaCB14, sig_cb14] = CB_2014_nga(M, T, Rrup, Rjb, Rx, W, Ztor, Zbot, delta, lambda, HW, VS30, Z25, Zhyp, region);
[Sa_bssa14, sig_bssa14] = BSSA_2014_nga(M, T, Rjb, lambda, region, Z1_0, VS30);
[SaPh20, ~, sig_Ph20,~] = Phung_2019h_NGAw2_TW(M, T, Rrup, Rjb, Rx, Ztor, delta, HW, lambda, VS30, Z1_0, region, fas);
 
FRO = 1; FSS = 0; FNO = 0; Finter = 0; Fintra = 0; Fas = 0;  Fma = 0;  FVS30 = 0;
[SaCh20,tau,phis2s,phiss] = Chao19H(T,M,Rrup,Ztor,VS30,Z1_0, FRO,FSS,FNO,Finter,Fintra,Fas,Fma,FVS30);
 
sigCh20 = sqrt(tau.^2 + phiss.^2);

median = [SaASK14, Sa_bssa14, SaCB14, SaCY14, SaCh20, SaPh20];
sigma = [sig_ask14, sig_bssa14, sig_cb14, sig_cy14, sig_Ph20, sigCh20];


fig=figure;
for k = 1:length(median)
    mu = median(k);
    sig = sigma(k);

    x = linspace(log(mu) - 4.0*sig, log(mu) + 4.0*sig, 100);

    % y = normpdf(x, log(SaASK14), sig_ask14);
    p = lognpdf(exp(x), mu, exp(sig));
    
    h(k) = semilogx(exp(x), p, 'LineWidth',2); hold on;
end
legend(h,{'ASK14','BSSA14','CB14','CY14','Ch20','Ph20'},'location','best');
title('M=7, R_{rup} = 30km, V_{S30}=760m/s');
xlabel('PGA, (g)'); grid on;
ylabel('Density');

fig_name = 'SampleGMMs.png';
saveas(fig, [dir_output, fig_name]);
 