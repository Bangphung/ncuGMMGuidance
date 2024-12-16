function [y, sigma] = Korea_GMM15(Tc, M, R)
    % Given only four periods

        T = [0, 0.01, 0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 2, 3];
    
%     if strcmp(model,'m0')
        a = [-3.16, -3.16, -3.15, -3.79,            -4.18,      -4.54, -4.73,    -4.95,         -5, -4.83, -4.84]; 
        b = [0.75, 0.75,    0.75, 0.85,                  0.94,   0.99,  1.02,    1.02,      0.98, 0.71, 0.57];
        c = [-0.72,-0.72,    -0.48, -0.45,               -0.52,   -0.64, -0.89,   -0.98,       -1.02, -0.94, -0.81]; 
        h = [3.7, 3.7,     3.9, 4.2,                       4.8,    5.9,  8.2,    8.9,      8.9, 7.8, 6.3]; 
        d = [-0.0034, -0.0034, -0.0039, -0.0035,             -0.0029,  -0.002, -0.0008,   -0.0005,    -0.0004, -0.0001, -0.0002];
        tau = [0.18, 0.18,  0.18, 0.17,              0.16,       0.15,  0.16,    0.16,         0.16, 0.16, 0.17]; 
        sig = [0.47, 0.47,  0.48, 0.43,              0.41,       0.39,  0.36,     0.34,        0.35, 0.34, 0.34];

%     elseif strcmp(model,'m1')
%         a = [-3.07,-3.07,	-4.27,	-5.15, -4.84];
%         b = [0.73, 0.73,	0.93,	0.95, 0.68];
%         c = [-0.76,	-0.76, -0.47,	-0.92, -0.91];
%         h = [1.7, 1.7,	2.3,	6.8, 6.4];
%         d = [-0.0029, -0.0029,	-0.0028,	-0.0003, 0.0002];
%         e = [0.326,	0.326, 2.77,	0.208, 0.21];
%         tau	= [0.17, 0.17,	0.14,	0.14, 0.15];
%         sig	= [0.34, 0.34,	0.31,	0.29, 0.29];
%         % sigT= [0.38,	0.34,	0.32];
%     end

    
    ip = find(T == Tc);
    logy = a(ip) + b(ip)*M + c(ip)*log(sqrt(R.^2 + h(ip)^2)) + d(ip)*R;% + e(ip)*s; 
    
    % unit of m/s^2;
    % convert to cm/s^2
    y = exp(logy)*100./980.6;
    
    sigma = sqrt(tau(ip)^2 + sig(ip)^2);

end