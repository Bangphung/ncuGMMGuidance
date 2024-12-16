clc; clear; close all;

Ti =  [0.01 0.02 0.03 0.04 0.05  0.075  0.1  0.12  0.15  0.17 0.2  0.25  0.3  0.4  0.5  0.75  1.0  1.5  2  3 4 5 7.5 10]';
nT = length(Ti);
M =  [4 5 6 7 8];
D = [1 1 10 10 30 30]; 
Ztor = 999; 
lambda = 0;
S = [270 760 270 760 270 760]; 
delta = 90;        
%% 
figure('position',[10 10 900 900]);
for j = 1:6 
     Rrup = D(j);
     Vs30 = S(j);
     
     subplot(3,2,j)
     
     Sav = arrayfun(@Phung_2019v_NGAw2_TW, [4 5 6 7 8].*ones(nT,1), repmat(Ti,1,5), Rrup*ones(nT,5), zeros(nT,5), zeros(nT,5),...
               Ztor*ones(nT,5),delta*ones(nT,5),zeros(nT,5),lambda*ones(nT,5),Vs30*ones(nT,5),...
               -999*ones(nT,5),ones(nT,5),zeros(nT,5),'UniformOutput',false);
           
     Sah = arrayfun(@Phung_2019h_NGAw2_TW, [4 5 6 7 8].*ones(nT,1), repmat(Ti,1,5), Rrup*ones(nT,5), zeros(nT,5), zeros(nT,5),...
               Ztor*ones(nT,5),delta*ones(nT,5),zeros(nT,5),lambda*ones(nT,5),Vs30*ones(nT,5),...
               -999*ones(nT,5),ones(nT,5),zeros(nT,5),'UniformOutput',false);
         
     pa = semilogx(Ti, cell2mat(Sav)./cell2mat(Sah),'-','LineWidth',1.8); hold on;
     
     axis([0.01 10  0  2.35]); 
     set(gca,'xtick',[0.01 0.1 1 3  10]);
     set(gca,'fontsize',12); grid on;
     set(gca,'Ytick',0:0.5:3);
     xlabel('Period, [sec]','fontsize',12);
     ylabel('V/H','fontsize',12);
     title(['SS; ','R_r_u_p=',num2str(Rrup),'km, ',...
            'V_S_3_0=',num2str(Vs30),'m/s']);
        if j ==1
            legend(pa,'M=4','M=5','M=6','M=7','M=8');
        end
end
%%

cm = [pa(1).Color; pa(2).Color; pa(3).Color; pa(4).Color;pa(5).Color];

%% 
figure('position',[10 10 900 900]);
for j = 1:6 
     Rrup = D(j);
     Vs30 = S(j);
     
     subplot(3,2,j)
     
     Sav = arrayfun(@Phung_2019v_NGAw2_TW, [4 5 6 7 8].*ones(nT,1), repmat(Ti,1,5), Rrup*ones(nT,5), zeros(nT,5), zeros(nT,5),...
               Ztor*ones(nT,5),delta*ones(nT,5),zeros(nT,5),lambda*ones(nT,5),Vs30*ones(nT,5),...
               -999*ones(nT,5),ones(nT,5),zeros(nT,5),'UniformOutput',false);
           
     Sah = arrayfun(@Phung_2019h_NGAw2_TW, [4 5 6 7 8].*ones(nT,1), repmat(Ti,1,5), Rrup*ones(nT,5), zeros(nT,5), zeros(nT,5),...
               Ztor*ones(nT,5),delta*ones(nT,5),zeros(nT,5),lambda*ones(nT,5),Vs30*ones(nT,5),...
               -999*ones(nT,5),ones(nT,5),zeros(nT,5),'UniformOutput',false);
         
     for k = 1:5
         pa(k)=loglog(Ti, cell2mat(Sav(:,k)),'-','color',cm(k,:),'LineWidth',1.8); hold on;
     
         pb=loglog(Ti, cell2mat(Sah(:,k)),'--','color',cm(k,:),'LineWidth',1.8); hold on;
     end
     
     
     axis([0.01 10  1e-4  3]); 
     set(gca,'xtick',[0.01 0.1 1 3  10]);
     set(gca,'fontsize',12); grid on;
     set(gca,'Ytick',[0.001 0.01 0.1 1]);
     xlabel('Period, [sec]','fontsize',12);
     ylabel('Sa, [g]','fontsize',12);
     title(['SS; ','R_r_u_p=',num2str(Rrup),'km, ',...
            'V_S_3_0=',num2str(Vs30),'m/s']);
        if j ==1
            legend(pa,'M=4','M=5','M=6','M=7','M=8','location','best');
        end
end