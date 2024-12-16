
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Variables
% M     = Moment Magnitude
% T     = Period (sec); 
%              
% Rrup   = Closest distance (km) to the ruptured plane
% Rjb   = Joyner-Boore distance (km); closest distance (km) to surface
%       projection of rupture plane
% Rx    =    Horizontal distance from top of rupture measured perpendicular 
%       to fault strike (km). See Figures a, b and c for illustation% Ztor  = Depth(km) to the top of ruptured plane
% delta = Fault dip angle (in degree)
% lamda = Rake angle      (in degree)
% region        = 1 for Taiwan
%               = 2 for California
%               = 4 for Japan
%               = 3 for others 

% Z10            = Basin depth (m); depth from the groundsurface to the
%                   1km/s shear-wave horizon.
%               = -999 if unknown
%               = 'na' if unknow
% Vs30          = shear wave velocity averaged over top 30 m in m/s
%               = ref: 1130
% Ztor          = depth to top of rupture , [km]
%               =-999 if specify average Ztor
% Output Variables
% Sa: Median spectral acceleration prediction
% sigma: logarithmic standard deviation of spectral acceleration
%          prediction based on tau and phi model of CY14 with the change in
%          phi1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Sa, yrefij, sigma_T, sigma_T2] = Phung_2019c2_Sub_TW(M, T, Rrup, Ztor, Vs30, Z10, Eqt)
   period = [0 0.01 0.02	0.03  0.04	0.05  0.075  0.1  0.12  0.15  0.17  0.2  0.25 0.30 0.4  0.5 ...
             0.75  1.0   1.5   2  3  4  5];
%% 
if length (T) == 1 && T == 1000 % Compute Sa and sigma with pre-defined period
    Sa=zeros(1,length(period));
    yrefij=zeros(1,length(period));

    sigma_T = zeros(1, length(T));
     sigma_T2 = zeros(1, length(T));
    for ip=1:length(period)
        [Sa(ip),yrefij(ip),sigma_T(ip),sigma_T2(ip)] = Phung_2018_sub(M, ip, Rrup, Ztor, Vs30, Z10, Eqt);
    end
else                            % Compute Sa and sigma with user-defined period
    Sa=zeros(1, length(T));
    yrefij = zeros(1, length(T));
    sigma_T = zeros(1, length(T));
     sigma_T2 = zeros(1, length(T));
    for i=1:length(T)
        Ti = T(i);
        if ( isempty(find(abs(period-Ti) < 0.0001))) % The user defined period requires interpolation
            T_low = max(period(period < Ti));
            T_high = min(period(period > Ti));
            ip_low  = find(period==T_low);
            ip_high = find(period==T_high);
            
            [Sa_low,yrefij_low,sig_low,sig_low2] = Phung_2018_sub(M, ip_low, Rrup, Ztor, Vs30, Z10, Eqt);
            [Sa_high,yrefij_high,sig_high,sig_high2] = Phung_2018_sub(M, ip_high, Rrup, Ztor, Vs30, Z10, Eqt);
           
            x = [log(T_low) log(T_high)];
            Y_sa = [log(Sa_low) log(Sa_high)];
            Y_yrefij = [yrefij_low yrefij_high];
            Y_sig = [sig_low sig_high];
             Y_sig2 = [sig_low2 sig_high2];
            
            Sa(i) = exp(interp1(x, Y_sa, log(Ti)));
            yrefij(i) = exp(interp1(x, Y_yrefij, log(T(i))));
            sigma_T(i) = interp1(x, Y_sig, log(T(i)));
            sigma_T2(i) = interp1(x, Y_sig2, log(T(i)));
        else
            ip_T = find(abs((period- Ti)) < 0.0001);
            [Sa(i),yrefij(i), sigma_T(i),sigma_T2(i)] = Phung_2018_sub(M, ip_T, Rrup, Ztor, Vs30, Z10, Eqt);
        end
    end
end


function [Sa, yrefij,Sig_SS,Sig_T] = Phung_2018_sub(M, ip, R_RUP, Ztor, Vs30, Z10,  Eqt)
%% period-indipendent fixing coefs
c4 = -2.1;
c4_a = -0.5;
c_RB = 50;
c2 = 1.06; 
%%
c1 = [-1.5659	-1.5131	-1.4626	-1.4145	-1.3689	-1.3256	-1.2280	-1.1452	-1.0892	-1.0217	-0.9884	-0.9636	-0.9497	-0.9308	-0.9218	-0.9384	-1.0733	-1.2193	-1.4437	-1.6454	-2.0106	-2.3871	-2.7714];
c1ss = [1.0392	1.0751	1.1075	1.1364	1.1617	1.1836	1.2231	1.2413	1.2413	1.2129	1.1762	1.1153	0.9926	0.8874	0.7548	0.6532	0.4641	0.2977	0.1023	0.0000	0.0000	0.0000	0.0000];

c3	= [0.199981482	-2.10416539	0.735340972	0.942876448	0.939915024	0.929699655	1.232745829	1.284670679	-9.172457382	0.950243111	0.982106563	1.022187278	1.043886732	1.061508105	1.23513106	1.328565104	1.506701615	1.636160785	1.873083499	1.939547698	2.010503041	2.010106353	2.053111123];
c_n	= [12.14866822	12.14866822	12.24803407	12.53378414	12.99189704	13.65075358	15.71447541	16.77262242	16.77563032	16.18679785	15.84314399	15.01467129	12.69643148	10.44981091	6.802216656	4.41069375	3.4064	3.1612	2.8078	2.4631	2.2111	1.9668	1.6671];
c_m	= [5.2316	5.3031	5.3734	5.4424	5.5103	5.5770	5.7384	5.8924	6.0103	6.1785	6.2849	6.4358	6.6646	6.8651	7.1814	7.3845	7.3830	7.3246	7.2838	7.3166	7.4292	7.5016	7.5328];

c5	= [12.4153	12.4153	12.3141	12.4522	12.1825	11.6987	10.4875	9.8032	9.6205	9.7344	9.9560	10.3810	11.1234	11.7692	12.7007	13.3302	14.6856	16.3825	20.5074	24.4969	30.2604	33.2809	34.5123];
c6	= [0.39264	0.39264	0.394	0.39936	0.40296	0.40384	0.40384	0.40384	0.40384	0.4036	0.40288	0.40128	0.39768	0.39352	0.38456	0.37656	0.366	0.36176	0.36008	0.36	0.36	0.36	0.36];

c_g1a = [-0.006936568	-0.006905602	-0.006881859	-0.007060489	-0.007167218	-0.007195002	-0.006622034	-0.006337324	-0.0063	-0.0066	-0.00679086	-0.007287566	-0.007136339	-0.007211224	-0.00690331	-0.006583927	-0.00680122	-0.006654415	-0.006659173	-0.006121232	-0.005478712	-0.006195018	-0.00581104];
c_g1b = [-0.00981353	-0.00980151	-0.009769199	-0.010047371	-0.010063113	-0.009850645	-0.008864109	-0.008609342	-0.0085	-0.0088	-0.009119758	-0.009688923	-0.0095747	-0.009443837	-0.008972698	-0.008485866	-0.008276172	-0.007539966	-0.007599027	-0.007226374	-0.006865062	-0.007358742	-0.007189709];

phi1	= [-0.510745033	-0.510415026	-0.502941955	-0.491366306	-0.474484696	-0.459984157	-0.446396645	-0.476282069	-0.4931516	-0.517925624	-0.532965478	-0.547665313	-0.565549294	-0.606451856	-0.653316566	-0.674933921	-0.796961941	-0.884871551	-0.958271065	-0.968084348	-0.96759396	-0.964753341	-0.923270348	-0.85471647	-0.770092758];
phi2	= [-0.1417	-0.1417	-0.1364	-0.1403	-0.1591	-0.1862	-0.2538	-0.2943	-0.3077	-0.3113	-0.3062	-0.2927	-0.2662	-0.2405	-0.1975	-0.1633	-0.1028	-0.0699	-0.0425	-0.0302	-0.0129	-0.0016	0.0000	0.0000	0.0000];
phi3	= [-0.007010 -0.007010	-0.007279	-0.007354	-0.006977	-0.006467	-0.005734	-0.005604	-0.005696	-0.005845	-0.005959	-0.006141	-0.006439	-0.006704	-0.007125	-0.007435	-0.008120	-0.008444	-0.007707	-0.004792	-0.001828	-0.001523	-0.001440	-0.001369	-0.001361];
phi4	= [0.102151	0.102151	0.108360	0.119888	0.133641	0.148927	0.190596	0.230662	0.253169	0.266468	0.265060	0.255253	0.231541	0.207277	0.165464	0.133828	0.085153	0.058595	0.031787	0.019716	0.009643	0.005379	0.003223	0.001134	0.000515];                                                                                                                      
%% Fixed Coefficients 
c7	= [0.00803536	0.00803536	0.007592927	0.007250488	0.007006048	0.006860143	0.007007726	0.007246641	0.007455965	0.00770271	0.007798775	0.007823121	0.00807121	0.008395901	0.00927498	0.010165826	0.012793392	0.013761922	0.013899933	0.012559337	0.009183764	0.004796976	0.001067909	-0.004234005	-0.006203311];
c7si = [0.0546	0.0861	0.0446	0.0409	0.0422	0.0440	0.0430	0.0494	-0.0046	0.0460	0.0502	0.0496	0.0517	0.0501	0.0456	0.0476	0.0429	0.0286	0.0188	0.0063	-0.0116	-0.0093	-0.0198];
c7ss = [-0.0126	-0.0407	-0.0122	-0.0089	-0.0077	-0.0069	0.0007	0.0014	-0.2351	-0.0132	-0.0131	-0.0124	-0.0132	-0.0159	-0.0125	-0.0119	-0.0166	-0.0210	-0.0174	-0.0182	-0.0152	-0.0077	-0.0001];

% c7si = [-0.5816	-1.0140	-0.6150	-0.5505	-0.5341	-0.5177	-0.4193	-0.4842	-5.1626	-0.7509	-0.7412	-0.6723	-0.6707	-0.6992	-0.5239	-0.5109	-0.3778	-0.3410	-0.2254	-0.2056	-0.0522	-0.0867	0.0010];
% c7si = [1.0245	1.2860	0.9727	0.9395	0.9537	0.9652	0.9927	1.0848	1.0667	1.1030	1.1238	1.0904	1.1023	1.0615	0.9144	0.9019	0.7313	0.5806	0.4974	0.3288	0.1018	0.1318	0.0067];

%%
c_g2	= [-0.007127092	-0.007127092	-0.007248737	-0.007327856	-0.007361759	-0.007360913	-0.007051574	-0.005719182	-0.00436511	-0.002649555	-0.001999512	-0.001254506	-0.00075041	-0.000447155	-0.000247246	-0.000416797	-0.001131462	-0.001741492	-0.002427965	-0.002705545	-0.004107346	-0.005776395	-0.007747849	-0.009141588	-0.012633296];
c_g3	= [4.225634814	4.225634814	4.230341898	4.236182109	4.250188668	4.303122568	4.446126947	4.610835161	4.723496543	4.878140618	4.981707053	5.066410859	5.21986472	5.32821979	5.201761713	5.187931728	4.877209058	4.63975087	4.571203643	4.425116502	3.6219035	3.48626393	3.277906342	3.074948086	3.074948086];
c_HM = [3.0956	3.0956	3.0963	3.0974	3.0988	3.1011	3.1094	3.2381	3.3407	3.4300	3.4688	3.5146	3.5746	3.6232	3.6945	3.7401	3.7941	3.8144	3.8284	3.8330	3.8361	3.8369	3.8376	3.8380	3.8380];
%%
tau = [0.372991982	0.372455184	0.37404339	0.386799496	0.400724323	0.415407087	0.425118035	0.416469357	0.402380557	0.382176649	0.371903621	0.357414054	0.337728984	0.359262599	0.397614543	0.426900573	0.466977537	0.495441465	0.487074394	0.477953808	0.436531699	0.449129802	0.46456534	0.505675779	0.448423418];
phi_SS = [0.4397	0.4388	0.4391	0.4451	0.4516	0.4555	0.4558	0.4497	0.4429	0.4382	0.4385	0.4395	0.4433	0.4498	0.4590	0.4703	0.4707	0.4643	0.4568	0.4521	0.4524	0.4461	0.4420	0.4177	0.3926];
phi_S2S = [0.3149	0.3149	0.3148	0.3223	0.3347	0.3514	0.3845	0.3935	0.3897	0.3713	0.3632	0.3503	0.3343	0.3324	0.3299	0.3319	0.3384	0.3480	0.3697	0.3826	0.3974	0.3983	0.3985	0.3878	0.3717];

sig_SS = [0.5766	0.5756	0.5768	0.5897	0.6038	0.6165	0.6233	0.6129	0.5984	0.5814	0.5750	0.5665	0.5573	0.5757	0.6073	0.6352	0.6630	0.6790	0.6678	0.6579	0.6287	0.6330	0.6412	0.6559	0.5960];
sig_TOT = [0.672899987	0.671374421	0.672061617	0.687372811	0.706122677	0.725510488	0.748191261	0.744953743	0.730667128	0.706545583	0.696399878	0.682643889	0.666628056	0.6810879	0.707107965	0.732357336	0.758959277	0.776876787	0.778458113	0.774272629	0.742033865	0.747199151	0.755236371	0.769755177	0.728570986];   
%% fmag
term6 = c2*(M-6);
term7 = (c2-c3(ip))/c_n(ip)*log(1+exp(c_n(ip)*(c_m(ip)-M)));
%% Distance Scaling and attenuation term
% del_c5 = dp(ip)*max(Ztor/50-20/50,0)*(Ztor>20)*(M<7.0); 
c_g1 = c_g1a(ip)*(Eqt==0) + c_g1b(ip)*(Eqt==1);
CNS = c5(ip)*cosh(c6(ip)*max(M-c_HM(ip),0));  
term8 =  c4*log(R_RUP + CNS);
term9 =  (c4_a-c4)*log(sqrt(R_RUP^2+c_RB^2));
term10 = (c_g1 + c_g2(ip)/(cosh(max(M-c_g3(ip),0))))*R_RUP;
%% Ztor term
if Eqt==1
   E_Ztor = (max(3.5384-2.60*max(M-5.8530,0),0)).^2;
else
   E_Ztor = 1.5*(max(2.7482-1.7639*max(M-5.5210,0),0)).^2;
end 
       
if Ztor == 999
   Ztor = E_Ztor;
   delta_ZTOR = 0;
else
   delta_ZTOR = Ztor - E_Ztor;
end
c7_b = c7si(ip)*(Eqt==0) + c7ss(ip)*(Eqt==1);
term4 = (c7(ip) + c7_b/cosh(2*max(M-4.5,0)))*delta_ZTOR;
%% main shock
   term1 = c1(ip) + c1ss(ip)*(Eqt==1);
   ln_yrefij = term1 + 0 + term6 + term7 + term8 + term9 + term10; 
   yrefij = exp(ln_yrefij);
%% Site response
term14 = phi1(ip)*min(log(Vs30/1130),0);
term15 = phi2(ip)*(exp(phi3(ip)*(min(Vs30,1130)-360))-exp(phi3(ip)*(1130-360)))*log((yrefij+phi4(ip))/phi4(ip));
%% Basin Depth   
%% median prediction
Sa = yrefij*exp(term14 + term15);
%% Variance Model-2 for Taiwan
Sig_SS = sqrt(phi_SS(ip)^2 + phi_S2S(ip)^2);
Sig_T = sqrt(tau(ip)^2 + phi_SS(ip)^2 + phi_S2S(ip)^2);
