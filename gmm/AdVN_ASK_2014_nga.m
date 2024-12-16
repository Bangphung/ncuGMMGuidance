% coded by Yue Hua, 5/19/10
%               Stanford University
%               yuehua@stanford.edu
% 
% Summary of the ASK14 Ground-Motion Relation for Active Crustal Regions
% Norman Abrahamson, Walter Silva and Ronnie Kamai
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Variables
% M     = Moment Magnitude
% T     = Period (sec); Use Period = -1 for PGV computation
%                 Use 1000 for output the array of Sa with original period
%                 (no interpolation)
% CRjb  = Centroid CRjb, assumed to be 999.9 here -> assume no aftershock
% Rrup   = Closest distance (km) to the ruptured plane
% Rjb   = Joyner-Boore distance (km); closest distance (km) to surface
%       projection of rupture plane
% Rx    =  Site coordinate (km) measured perpendicular to the fault strike
%       from the fault line with down-dip direction to be positive
% Ry0   = Horizontal distance off the end of the rupture measured oarallel
%       to strike
% Ztor  = Depth(km) to the top of ruptured plane
% delta = Fault dip angle (in degree)
% lamda = Rake angle      (in degree)
% fas   = Flag for aftershocks
% HW    = Falg for hanging wall sites
% W     = Down-dip rupture width (km) 
% region        = 0 for global
%               = 1 for California
%               = 2 for Japan
%               = 3 for China 
%               = 4 for Italy 
%               = 5 for Turkey
%               = 6 for Taiwan
% Z10            = Basin depth (km); depth from the groundsurface to the
%                   1km/s shear-wave horizon.
%               = 999 if unknown
%               = 'na' if unknow
% Vs30          = shear wave velocity averaged over top 30 m in m/s
%               = ref: 1130
% FVS30         = 1 for Vs30 is inferred from geology
%               = 0 for measured  Vs30

% Output Variables
% Sa: Median spectral acceleration prediction
% sigma: logarithmic standard deviation of spectral acceleration
%          prediction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Sa, sigma] = AdVN_ASK_2014_nga(M, T, Rrup, Rjb, Rx, Ry0, Ztor, delta, lambda, fas, HW, W, Z10, Vs30, SCL, region)

period = [0.01	0.02	0.03	0.05	0.075	0.1	0.15	0.2	0.25	0.3	0.4	0.5	0.75	1	1.5	2	3	4	5	6	7.5	10	0	-1];
frv = lambda >= 30 & lambda <= 150; % frv: 1 for lambda between 30 and 150, 0 otherwise
fnm = lambda >= -150 & lambda <= -30; % fnm: 1 for lambda between -150 and -30, 0 otherwise

if Ztor == 999
    if frv == 1
        Ztor = max(2.704 - 1.226 * max(M-5.849,0),0)^2;
    else
        Ztor = max(2.673 - 1.136 * max(M-4.970,0),0)^2;
    end
end

if W == 999
    W = min(18/sin(degtorad(delta)),10^(-1.75+0.45*M));
end

if Z10 == 999
    if region == 2
        Z10 = exp(-5.23 / 2 * log((Vs30 ^ 2 + 412 ^ 2) / (1360 ^ 2 + 412 ^ 2))) / 1000;
    else %'non-japanese
        Z10 = exp((-7.67 / 4) * log((Vs30 ^ 4 + 610 ^ 4) / (1360 ^ 4 + 610 ^ 4))) / 1000;
    end
end

if length (T) == 1 && T == 1000; % Compute Sa and sigma with pre-defined period
    Sa=zeros(1,length(period)-2);
    sigma=zeros(1,length(period)-2);
    period1=period(1:end-2);
    for ip=1:length(period)-2
        [Sa(ip),sigma(ip)]=ASK_2014_sub_1(M, ip, Rrup, Rjb, Rx, Ry0, Ztor, delta, frv, fnm, fas, HW, W, Z10, Vs30, SCL, region);
    end
else                            % Compute Sa and sigma with user-defined period
    Sa=zeros(1, length(T));
    sigma=zeros(1, length(T));
    period1=T;
    for i=1:length(T)
        Ti = T(i);
        if ( isempty(find(abs(period-Ti) < 0.0001))) % The user defined period requires interpolation
            T_low = max(period(period < Ti));
            T_high = min(period(period > Ti));
            ip_low  = find(period==T_low);
            ip_high = find(period==T_high);
            
            [Sa_low, sigma_low] = ASK_2014_sub_1(M, ip_low, Rrup, Rjb, Rx, Ry0, Ztor, delta, frv, fnm, fas, HW, W, Z10, Vs30, SCL,  region);
            [Sa_high, sigma_high] = ASK_2014_sub_1(M, ip_high, Rrup, Rjb, Rx, Ry0, Ztor, delta, frv, fnm, fas, HW, W, Z10, Vs30, SCL,  region);
            x = [log(T_low) log(T_high)];
            Y_sa = [log(Sa_low) log(Sa_high)];
            Y_sigma = [sigma_low sigma_high];
            Sa(i) = exp(interp1(x, Y_sa, log(Ti)));
            sigma(i) = interp1(x, Y_sigma, log(Ti));
        else
            ip_T = find(abs((period- Ti)) < 0.0001);
            [Sa(i), sigma(i)] = ASK_2014_sub_1(M, ip_T, Rrup, Rjb, Rx, Ry0,Ztor, delta, frv, fnm, fas, HW, W, Z10, Vs30, SCL,  region);
        end
    end
end

function [Sa, sigma]=ASK_2014_sub_1 (M, ip, R_RUP, R_JB, Rx, Ry0, Ztor, delta, F_RV, F_NM, F_AS, HW, W, Z10, Vs30, SCL, region)

%% Coefficients
T = [0.01	0.02	0.03	0.05	0.075	0.1	0.15	0.2	0.25	0.3	0.4	0.5	0.75	1	1.5	2	3	4	5	6	7.5	10	0	-1];
Vlin =[	660.0000	680.0000	770.0000	915.0000	960.0000	910.0000	740.0000	590.0000	495.0000	430.0000	360.0000	340.0000	330.0000	330.0000	330.0000	330.0000	330.0000	330.0000	330.0000	330.0000	330.0000	330.0000	660.0000	330.0000];
b	=[-1.4700	-1.4590	-1.3900	-1.2190	-1.1520	-1.2300	-1.5870	-2.0120	-2.4110	-2.7570	-3.2780	-3.5990	-3.8000	-3.5000	-2.4000	-1.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	-1.4700	-2.0200];
n	=[1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000];
M1	=[6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.8200	6.9200	7.0000	7.0600	7.1450	7.2500	6.7500	6.7500];
c	=[2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2400.0000];
c4	=[4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000];



% a1VN = 0.*[-0.9007	-0.9248	-0.7549	-0.6998	-0.8249	-0.7966	-0.8222	-1.0377	-1.2173	-1.4279	-1.6888	-1.8215	-1.7990	-1.7925	-1.9903	-1.8560];
% a6VN = 0.*[-0.9178	-0.9883	-1.0352	-1.0588	-0.9831	-0.9180	-0.9271	-0.9441	-0.9409	-0.9280	-0.9125	-0.8740	-0.8018	-0.7743	-0.7960	-0.7947];


% a1VN = 0.*[-1.1362	-1.1010	-0.9395	-0.8211	-0.9864	-1.1360	-1.2002	-1.3422	-1.4687	-1.5667	-1.7715	-1.9280	-1.8944	-1.7718	-1.8958	-1.7597];

% a1VN = 0.*[-1.1620	-1.1042	-0.8900	-0.8071	-0.9920	-1.1211	-1.1907	-1.3505	-1.4676	-1.6281	-1.7836	-1.9373	-1.8896	-1.7740	-1.8811	-1.7668];


a1VN = 1.*[-0.6015	-0.5410	-0.3243	-0.2367	-0.4189	-0.5463	-0.6124	-0.7742	-0.8970	-1.0673	-1.2549	-1.4316	-1.4240	-1.3335	-1.4473	-1.3330];

a6VN = 0.5*[-0.9433	-0.9473	-0.9510	-0.9576	-0.9643	-0.9693	-0.9742	-0.9723	-0.9636	-0.9480	-0.8937	-0.8550	-0.7869	-0.7448	-0.7448	-0.7448];

a7VN = 0.*[0.4151	0.4226	0.4335	0.4204	0.4284	0.4491	0.4489	0.4332	0.4348	0.4064	0.4051	0.3926	0.3744	0.3622	0.3880	0.3928];

a4VN = 1.*[-0.3494	-0.3475	-0.3456	-0.3426	-0.3396	-0.3376	-0.3361	-0.3387	-0.3466	-0.3522	-0.3525	-0.3434	-0.3004	-0.2665	-0.2232	-0.1965];

a17VN = 1.*[-0.00680	-0.00707	-0.00731	-0.00769	-0.00799	-0.00811	-0.00780	-0.00688	-0.00576	-0.00467	-0.00302	-0.00181	-0.00044	0.00000	0.00000	-0.00032];


a1	=[0.5870	0.5980	0.6020	0.7070	0.9730	1.1690	1.4420	1.6370	1.7010	1.7120	1.6620	1.5710	1.2990	1.0430	0.6650	0.3290	-0.0600	-0.2990	-0.5620	-0.8750	-1.3030	-1.9280	0.5870	5.9750];
a2	=[-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7650	-0.7110	-0.6340	-0.5290	-0.7900	-0.9190];
a3	=[0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750];
a4	=[-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000];
a5	=[-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100];
a6	=[2.1541	2.1461	2.1566	2.0845	2.0285	2.0408	2.1208	2.2241	2.3124	2.3383	2.4688	2.5586	2.6821	2.7630	2.8355	2.8973	2.9061	2.8888	2.8984	2.8955	2.8700	2.8431	2.1541	2.3657];
a7	=[0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000];
a8	=[-0.0150	-0.0150	-0.0150	-0.0150	-0.0150	-0.0150	-0.0220	-0.0300	-0.0380	-0.0450	-0.0550	-0.0650	-0.0950	-0.1100	-0.1240	-0.1380	-0.1720	-0.1970	-0.2180	-0.2350	-0.2550	-0.2850	-0.0150	-0.0940];
a10	=[1.7350	1.7180	1.6150	1.3580	1.2580	1.3100	1.6600	2.2200	2.7700	3.2500	3.9900	4.4500	4.7500	4.3000	2.6000	0.5500	-0.9500	-0.9500	-0.9300	-0.9100	-0.8700	-0.8000	1.7350	2.3600];
a11	=[0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000];
a12	=[-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.2000	-0.2000	-0.2000	-0.1000	-0.1000];
a13	=[0.6000	0.6000	0.6000	0.6000	0.6000	0.6000	0.6000	0.6000	0.6000	0.6000	0.5800	0.5600	0.5300	0.5000	0.4200	0.3500	0.2000	0.0000	0.0000	0.0000	0.0000	0.0000	0.6000	0.2500];

M2 = 5.0;
%% Site Effects
Ps = [1.1367, 1.5067, 2.1887, 2.9436, 2.4314];
c1A	= [-0.3718	-0.3527	-0.3384	-0.3235	-0.3288	-0.3641	-0.4757	-0.4814	-0.4212	-0.3751	-0.3530	-0.3652	-0.3889	-0.3507	-0.1897	-0.1321];
c1B	= [0.0780	0.0533	0.0351	0.0173	0.0286	0.0787	0.2371	0.2382	0.1578	0.1295	0.1632	0.1644	0.3281	0.3046	0.0405	-0.0763];
c2	= [0.6376	0.6404	0.6456	0.6622	0.6927	0.7323	0.8457	0.7572	0.5684	0.4142	0.2904	0.2462	0.2347	0.1733	0.1944	0.2379];
c3	= [0.8964	0.8465	0.8065	0.7553	0.7417	0.7892	1.0158	1.1758	1.2049	1.1808	0.9968	0.6922	0.4299	0.3020	0.2754	0.3274];
c4	= [0.9379	0.9054	0.8793	0.8449	0.8312	0.8572	0.9603	1.0353	1.0324	1.0399	1.0963	1.1846	1.2943	1.2381	0.7534	0.5498];

if ~ismember(SCL,[1,2,3,4,5])
    %SCL = 1;
    fSCL = c1A(ip)*log(Ps(1));
else
    if SCL == 1
       fSCL = c1A(ip)*log(Ps(1));
    elseif SCL ==2
         fSCL = c1B(ip)*log(Ps(2));
    elseif SCL ==3
         fSCL = c2(ip)*log(Ps(3));
    elseif SCL ==4
         fSCL = c3(ip)*log(Ps(4));
    else 
         fSCL = c4(ip)*log(Ps(5));
    end
end

%% Term f1 - Basic form
if M > 5
    c4m = c4(ip);
elseif M > 4 && M <= 5
    c4m = c4(ip)-(c4(ip)-1)*(5-M);
else
    c4m = 1;
end

R = sqrt(R_RUP^2+ c4m^2);

if M > M1(ip)
    f1 = a1(ip) + a1VN(ip) + a5(ip)*(M - M1(ip)) + a8(ip)*(8.5 - M)^2 + (a2(ip) + a3(ip)*(M - M1(ip)))*log(R) + a17VN(ip)*R_RUP;
elseif M > M2 && M <= M1(ip)
    f1 = a1(ip) + a1VN(ip) + a4(ip)*(M - M1(ip)) + a8(ip)*(8.5 - M)^2 + (a2(ip) + a3(ip)*(M - M1(ip)))*log(R) + a17VN(ip)*R_RUP;
else
    f1 = a1(ip) + a1VN(ip) + a4(ip)*(M2 - M1(ip)) + a8(ip)*(8.5 - M2)^2 + (a6(ip) + a6VN(ip))*(M - M2) + (a7(ip) + a7VN(ip))*(M - M2)^2 + ...
        (a2(ip) + a3(ip)*(M2 - M1(ip)))*log(R) + a17VN(ip)*R_RUP;
end

%% term f4 - Hanging wall model
R1 = W * cos(degtorad(delta));
R2 = 3 * R1;
Ry1 = Rx * tan(degtorad(20));
h1 = 0.25;
h2 = 1.5;
h3 = -0.75;

if delta > 30
    T1 = (90- delta)/45;
else
    T1 = 60/45;
end

a2hw = 0.2;

if M > 6.5
    T2 = 1 + a2hw * (M - 6.5);
elseif M > 5.5
    T2 = 1 + a2hw * (M - 6.5) - (1 - a2hw) * (M - 6.5)^2;
else
    T2 = 0;
end

if Rx <= R1
    T3 = h1 + h2*(Rx/R1) + h3*(Rx/R1)^2;
elseif Rx < R2
    T3 = 1 - (Rx - R1)/(R2 - R1);
else 
    T3 = 0;
end

if Ztor < 10
    T4 = 1 - Ztor^2/100;
else
    T4 = 0;
end

if Ry0 == 999 || Ry0 == 0
    if R_JB == 0
        T5 = 1;
    elseif R_JB <30
        T5 = 1 - R_JB/30;
    else
        T5 = 0;
    end
else
    if Ry0 - Ry1 <= 0
        T5 = 1;
    elseif Ry0 - Ry1 < 5
        T5 = 1- (Ry0-Ry1)/5;
    else
        T5 = 0;
    end
end

if HW == 1
    f4 = a13(ip) * T1 * T2 * T3 * T4 * T5;
else
    f4 = 0;
end

%% Term f6 - Depth to top rupture model
% if Ztor < 20
%     f6 = a15(ip)*Ztor/20;
% else
%     f6 = a15(ip);
% end
f6 = a4VN(ip)*min(Ztor/20,1);
%% Term: f7 and f8 - Style of Faulting

if M > 5
    f7 = a11(ip);
    f8 = a12(ip);
elseif M >= 4 && M <= 5
    f7 = a11(ip)*(M - 4);
    f8 = a12(ip)*(M - 4);
else 
    f7 = 0;
    f8 = 0;
end

if T(ip) <= 0.5
    V1 = 1500;
elseif T(ip) < 3
    V1 = exp(-0.35*log(T(ip)/0.5)+log(1500));
else 
    V1 = 800;
end  

if Vs30 < V1
    Vs30s = Vs30;
else
    Vs30s = V1;
end

if 1180 >= V1
    Vs30star1180 = V1;
else
    Vs30star1180 = 1180;
end

%% term  Regional:

%% Term f5 - site response model

%Sa 1180
f5_1180 = (a10(ip) + b(ip) * n(ip)) * log(Vs30star1180 / Vlin(ip));

Sa1180 = exp(f1 + f6 + F_RV*f7 + F_NM * f8 +  HW * f4 + f5_1180);

if Vs30 >= Vlin(ip)
    f5 = (a10(ip)+ b(ip)*n(ip))*log(Vs30s/Vlin(ip));
else 
    f5 = a10(ip)*log(Vs30s/Vlin(ip)) - b(ip)*log(Sa1180 + c(ip)) + b(ip)*log(Sa1180 + c(ip)*(Vs30s/Vlin(ip))^n(ip));
end

%% Term f10 - soil depth model

%% Sa

lnSa = f1 + 1*f6 + F_RV*f7 + F_NM * f8  +  fSCL;
 
Sa= exp(lnSa);
 
 
%% Standard deviation 

phi = 0;

tau = 0;

sigma = sqrt(phi^2+ tau^2);



    













    







