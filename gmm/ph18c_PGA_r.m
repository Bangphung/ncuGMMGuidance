function PGA_r = ph18c_PGA_r(M, Rrup, Vs30, Ztor, eqt, flag)
% M = xdata(:,1); Rrup = xdata(:,2); Vs30 = xdata(:,3); 
% Ztor = xdata(:,4); reg = xdata(:,5); eqt = xdata(:,6);
%% Fixing Coefs
c4 = 10;
a3 = 0.1; 
a9 = 0.25;
Tc = 0;
n = 1.18;
c = 1.88; 
%% free coefs  

period = [0 0.01 0.02  0.05  0.075  0.1  0.15  0.20  0.25  0.3  0.4  0.5  0.6  0.75  ... 
         1  1.5  2  2.5  3  4  5];    
Vlin	= [865.1	865.1	865.1	1053.5	1085.7	1032.5	877.6	748.2	654.3	587.1	503	456.6	430.3	410.5	400	400	400	400	400	400	400];
b	= [-1.186	-1.186	-1.186	-1.346	-1.471	-1.624	-1.931	-2.188	-2.381	-2.518	-2.657	-2.669	-2.599	-2.401	-1.955	-1.025	-0.299	0	0	0	0];
%%
a5 = [0.03849929	0.04033665	0.04190178	0.04509359	0.04623140	0.04819708	0.04325090	0.03692059	0.06597319	0.06197944	0.06979644	0.08783791	0.09612877	0.10612877	0.22744484	0.16136621	0.22767232	0.27153377	0.28822087	0.32589322	0.30383949];
a13 = [-0.0256568	-0.0259617	-0.0262528	-0.0270426	-0.0276048	-0.0280794	-0.0287650	-0.0291017	-0.0290970	-0.0287552	-0.0269993	-0.0235859	-0.0180673	-0.0150673	-0.0031849	-0.0031849	-0.0031849	-0.0031849	-0.0031849	-0.0031849	-0.0031849];
Mref = [7.68	7.68	7.68	7.71	7.77	7.77	7.78	7.72	7.62	7.54	7.42	7.38	7.36	7.32	7.25	7.25	7.25	7.25	7.25	7.25	7.25];   
a2  = [-1.552846733	-1.554174269	-1.555152194	-1.556049687	-1.554562252	-1.551165488	-1.539140832	-1.520764226	-1.489051706	-1.464118878	-1.414761429	-1.383170353	-1.360022278	-1.313716982	-1.236841977	-1.100570482	-0.990254902	-0.896093506	-0.818199517	-0.730697376	-0.734817372];
a14 = [-0.011876681	-0.012409284	-0.016872732	-0.08510905	-0.118005772	-0.171218187	-0.124720279	-0.120958201	-0.116255248	-0.077408811	-0.054966213	-0.034173086	-0.06069315	-0.039053473	0.017806808	-0.005705423	0.053155037	0.068765677	0.071577687	0.042405486	0.054712361];
%%
del_a1 = [1.141899742	1.152006702	1.154339185	1.515303155	1.904431142	1.945456526	1.787100626	1.562515125	1.356740101	1.206013896	0.760110718	0.431629072	0.214072689	-0.01782956	-0.204991951	-0.382378342	-0.352611734	-0.228719047	-0.16534756	0.010400184	0.135306871];
del_a4	=  [0.328613796	0.352192886	0.367677006	0.452541112	0.513739193	0.499522237	0.45427803	0.363869484	0.314270529	0.282854636	0.176621233	0.077343234	0.044354701	-0.012346815	-0.083175191	-0.226959428	-0.206168741	-0.163340854	-0.154799226	-0.112793712	-0.003105582];

a6_jp = [-0.006794362	-0.006817094	-0.006816137	-0.007285287	-0.007702243	-0.007674043	-0.00782682	-0.007547403	-0.007323965	-0.006976346	-0.006143907	-0.005504091	-0.004739568	-0.004285202	-0.003957479	-0.003379651	-0.003469497	-0.003409644	-0.003614492	-0.003749858	-0.003243671];
%a12_jp = [-0.7516020	-0.7500948	-0.7307185	-0.4831132	-0.3413025	-0.4948081	-0.8669192	-1.0634892	-1.1789740	-1.2253631	-1.2073943	-1.1299835	-1.0859780	-1.0233555	-0.9766258	-0.9437327	-0.8880212	-0.8546554	-0.7803988	-0.6937169	-0.6499105];
%a8_jp = [0.002650526	0.002810000	0.002376243	0.000548389	0.003562477	0.006709533	0.007165476	0.004341874	0.004771685	0.005747192	0.005294191	0.002624524	0.00028237	-0.001824156	-0.003649632	-0.005143864	-0.005623872	-0.004503918	-0.005225851	-0.006465579	-0.004203115];
%%
a1	= [4.46422	4.48194	4.50464	4.62523	4.79044	4.92617	5.09496	5.18110	5.13520	5.07746	4.91507	4.79323	4.66100	4.40139	3.74331	2.82226	1.88325	1.01478	0.41174	-0.51990	-0.92355];

a4	= [0.441987425	0.442328411	0.436082809	0.363261281	0.319469336	0.325896968	0.350561168	0.401101385	0.440779304	0.486141867	0.593885499	0.719248494	0.848115514	0.96522535	1.174894012	1.360979471	1.38307024	1.382803719	1.391730855	1.36799257	1.379913313];
a7 = [0.681875399	0.679916748	0.697566565	1.036117503	1.229932174	1.534436374	1.263659339	1.177341019	1.046187707	0.783076472	0.56945524	0.437054838	0.497762038	0.262897537	-0.126754198	-0.121313834	-0.496676074	-0.513466165	-0.600511124	-0.425097306	-0.545699974+0.0171];
a6tw = [-0.000639314	-0.000607826	-0.000577165	-0.000490096	-0.00042309	-0.000361045	-0.000251542	-0.000161014	-8.90E-05	-3.54E-05	1.45E-05	-4.21E-05	-7.26E-05	-0.000119483	-0.000199116	-0.000362196	-0.000611716	-0.000869846	-0.00106641	-0.001185004	-0.0009885];

Si12tw = [0.99033	0.99037	0.99282	1.31918	1.50198	1.63759	1.89340	2.08726	2.23477	2.34640	2.48809	2.50062	2.40243	2.07480	1.48948	0.38514	-0.41535	-0.75365	-0.73641	-0.69136	-0.70127];
Si12jp = [0.94642	0.94648	0.95388	1.31215	1.49926	1.61904	1.77120	1.87989	1.95713	2.03887	2.24946	2.29652	2.23798	2.05421	1.55696	0.57387	-0.26271	-0.63003	-0.54192	-0.49677	-0.42910];

a10	= [0.016025291	0.017193978	0.01828222	0.020842842	0.022162011	0.022757257	0.021797461	0.020180594	0.018556649	0.016978648	0.014555899	0.012627818	0.011191399	0.009211979	0.006851124	0.003814084	0.001733925	0	0	0	0];
a11 = [0.014951807	0.014930723	0.01491224	0.01487185	0.014854753	0.014852358	0.014893477	0.015004213	0.015194298	0.01540766	0.015952307	0.016437613	0.01652538	0.016212382	0.015784785	0.01399451	0.011927777	0.009749305	0.007785629	0.00494863	0.003408571];
%%
it = find (period==Tc);
%%
if flag ==0
   a1 = a1(it) + del_a1(it);
   a4 = a4(it) + del_a4(it);
   a6 = a6_jp(it);
   a12 = Si12jp(1);
else
   a1 = a1(it); 
   a4 = a4(it);
   a6 = a6tw(it);
   a12 = Si12tw(1);
end
%%
for i = 1:length(M)  
%% Magnitude Scaling
if M(i) <= Mref(it)
   fmag(i,1) =  a4*(M(i) - Mref(it)) + a13(it)*(10-M(i))^2;
else
   fmag(i,1) =  a5(it)*(M(i) - Mref(it)) + a13(it)*(10-M(i))^2;
end  
%% FZtor
if eqt ==0 
   f_ztor(i,1) = a10(it)*(min(Ztor(i),40)-20);
else
   f_ztor(i,1) = a11(it)*(min(Ztor(i),80)-40); 
end
%% Path Scaling
X(i,1) = a7(it)*(eqt(i)==1) + (a2(it) + a14(it)*(eqt(i)==1) + a3*(M(i) -7.8))*log(Rrup(i) + c4*exp(a9*(M(i)-6)))+ a6*Rrup(i);
%% The site function:
Vs30 = min(Vs30,1000);
fsite(i,1) =  a12*log(Vs30/Vlin(1))+ b(1)*n*log(Vs30/Vlin(1));
%%
%% 
PGA_r(i,1) = exp(a1 + fmag(i,1) + X(i,1) +  f_ztor(i,1) + fsite(i,1));
end
