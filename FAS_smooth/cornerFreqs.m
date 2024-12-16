function [hp_freq, lp_freq, fas_sig, loc, SNR] = cornerFreqs(x, dt, varargin)
%   AUTODETECT BANDPASS FILTER CORNER FREQUENCIES
%
%   cornerFreqs automatically detects appropriate bandpass filter corner
%   frequencies to be used for Butterworth filtering by using the seismic
%   signal's frequency content.
%
%   MOTIVATION:
%
%   Processing of seismic waveforms often requires bandpass filtering.
%   Selection of filter corner frequencies has been not only a manual
%   process but also subjective. There is a need for automatically
%   detecting corner frequencies for processing a large number of seismic
%   recordings.
%
%   ALGORITHM:
%
%   First, "PphasePicker" function (Kalkan, 2016) is used to determine
%   p-phase arrival time (event onset) to get background noise. Next,
%   Fourier amplitude spectra for noise and signal are calculated.
%   These two spectra are smoothed using "smoothSpectra"
%   function. Finally, intersection points of the smoothed spectra within
%   low-pass and high-pass frequency regions are searched to determine the
%   appropriate corner frequencies to be used for bandpass filtering.
%
%   High-pass region is defined between 0.1 Hz and 1 Hz. If no intersection
%   point detected, default value of 0.1 Hz is used.
%
%   Low-pass region is defined between the characteristic frequency
%   of the recording instrument (fc) (often 25 Hz) and Nyquist (half of
%   sampling frequency of the waveform data). If no intersection point
%   detected, 80% of Nyquist is used as the low-pass corner.
%
%   This code uses the following external functions:
%
%   [1] PphasePicker.m --> This function computes P-Phase onset time,
%       also available at MatLAB FEX
%
%   [2] smoothSpectra.m --> This function smooth FAS using Konno-Ohmachi
%       window, also available at MatLAB FEX
%
%   USAGE:
%
%   [hp_freq, lp_freq] = cornerFreqs(x,dt)
%
%   STATIC INPUT:
%
%            x = broadband velocity or acceleration data in
%                single-column format (1xn) or (nx1)
%           dt = sampling interval in second (e.g., 0.005)
%
%   VALID PROP_NAME / PROP_VAL PAIRS:
%   -----------------------------------------
%   'plot_name'    --> [text]-[default: None]
%   'plot_path'    --> [text]-[default: None]
%   'debug'        --> [text]-[default: False]
%
%   OUTPUT:
%
%      hp_freq = high-pass corner frequency in Hz
%      lp_freq = low-pass  corner frequency in Hz
%
%   EXAMPLES:
%
%   see demo.m file
%
%   REQUIREMENTS:
%
%   cornerFreqs function does not require any MatLAB toolbox.
%
%   ACKNOWLEDGEMENT:
%
%   In preparing this function, I benefitted from Curve Intersections
%   (InterX.m) function written by NS, which is available at MathWorks FEX.
%
%   REFERENCE:
%
%   Kalkan, E. (2016). ?An Automatic P-phase Arrival Time Picker,? Bulletin
%   of Seismological Society of America,106(3): 971-986, doi:
%   10.1785/0120150111.
%
%   THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED
%   WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
%   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
%   NO EVENT SHALL THE COPYRIGHT OWNER BE LIABLE FOR ANY DIRECT, INDIRECT,
%   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
%   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
%   OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
%   TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
%   USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
%   DAMAGE.
%
%   Written by Dr. Erol Kalkan, P.E. (kalkan76@gmail.com)
%   URL: www.erolkalkan.com
%   $Revision: 1.0.7 $  $Date: 2019/02/11 12:00:00 $
%
%% DEFAULT PROPERTIES
if (nargin >= 2)
    plot_name = '';
    plot_path = '';
    debug = 'False';
else
    error('cornerFreqs: First and second arguments must be a waveform and',...
          'sampling rate')
end

%% USER-DEFINED PROPERTIES
if (nargin > 2)
    v = varargin;
    nv = nargin-2;
    if ~rem(nv,2) == 0
        error(['cornerFreqs: Arguments after sampling rate must appear in ',...
            'property name/val pairs'])
    end
    for n = 1:2:nv-1
        name = lower(v{n});
        val = v{n+1};
        switch name
            case 'plot_name'
                plot_name = val;
            case 'plot_path'
                plot_path = val;
            case 'debug'
                debug = val;
            otherwise
                error('cornerFreqs: Property name not recognized')
        end
    end
end

%% Main
% enforce input as column vector
if isrow(x); x = x'; end

% run PphasePicker to get onset location
type = 'SM';    % no bandpass filtering
pflag = 'n';    % plot waveform and p-phase onset
Tn = 0.1;       % undamped natural period of oscillator in second
xi = 0.6;       % damping ratio
nbins = 100;    % histogram bin size
o = 'to_peak';  % 'to_peak' to take segment of waveform from beginning to peak

[loc, snr_db, ~] = PphasePicker(x, dt, type, pflag, Tn, xi, nbins, o, plot_path, plot_name);

% detrend signal by removing mean
signal = detrend(x, 'constant');

% extract pre-event segment, which we assume as background "noise"
noise = signal(1:round(loc/dt));
Dt = (length(noise) - 1)*dt; 

if Dt < 20
   npt = 30/dt + 1;
   noise = signal(end:-1:end-npt);
   disp(['len noise = ', num2str(length(noise))]);
end

% compute FAS and smooth FAS for noise
[fas_noi] = computeFAS(noise,dt);
[fas_noi.smooth] = smoothSpectra(fas_noi.norm,'w',50);

% compute FAS and smooth FAS for signal
[fas_sig] = computeFAS(signal,dt);
[fas_sig.smooth] = smoothSpectra(fas_sig.norm,'w',50);

% fc = characteristic frequency of most of the recording instruments
fc = 25;

% low-pass search: use fc as the lower limit of search for
% intersection points between two spectra
lim = fc;
% find index of freq arrays nearest to fc
[~, index_sig] = min(abs(fas_sig.freq-lim));
[~, index_noi] = min(abs(fas_noi.freq-lim));

% find index of freq arrays nearest to Nyquist
[~, index_sig_Ny] = min(abs(fas_sig.freq-0.5/dt));
[~, index_noi_Ny] = min(abs(fas_noi.freq-0.5/dt));

% get intersection points between signal and noise spectra (lp means
% low-pass)
    Plp = [];
% Plp = getInterX(fas_sig.freq(index_sig:index_sig_Ny),fas_sig.smooth(index_sig:index_sig_Ny), ...
%       fas_noi.freq(index_noi:index_noi_Ny),fas_noi.smooth(index_noi:index_noi_Ny));

 xq = sort(unique([fas_sig.freq(index_sig:index_sig_Ny), fas_noi.freq(index_noi:index_noi_Ny)]));
 fas_sig_smooth = interp1(fas_sig.freq(index_sig:index_sig_Ny),fas_sig.smooth(index_sig:index_sig_Ny),xq);
 fas_noi_smooth = interp1(fas_noi.freq(index_noi:index_noi_Ny),fas_noi.smooth(index_noi:index_noi_Ny),xq);
 [~, idmin] =  min(abs(fas_sig_smooth - fas_noi_smooth));

 Php = xq(idmin);  


if ~isempty(Plp)
    lp_freq = Plp(1,1);
    fprintf('cornerFreqs: low-pass frequency (Hz)......... %5.0f\n',lp_freq);
else % if no intersection
    fprintf('cornerFreqs: no low-pass corner frequency detected use %2.0f (Hz) (80 percent of Nyquist)\n', 0.4/dt);
    % lp_freq = 0.4/dt;
    lp_freq = 20;
end

% high-pass search: % use 1 Hz as the higher limit of search for
% intersection points between two spectra
lim = 1;
clear index_sig index_noi
[~, index_sig] = min(abs(fas_sig.freq-lim));
[~, index_noi] = min(abs(fas_noi.freq-lim));

hp_freq_min = 0.1; % lower limit for hp_freq search
% find index of freq arrays nearest to lower limit
[~, index_sig_min] = min(abs(fas_sig.freq-hp_freq_min));
[~, index_noi_min] = min(abs(fas_noi.freq-hp_freq_min));

% get intersection points between signal and noise spectra (hp means
% high-pass)
Php = [];
if index_noi > 1
%      Php = getInterX(fas_sig.freq(index_sig_min:index_sig), fas_sig.smooth(index_sig_min:index_sig), ...
%          fas_noi.freq(index_noi_min:index_noi),fas_noi.smooth(index_noi_min:index_noi));

     xq = sort(unique([fas_sig.freq, fas_noi.freq]));
     fas_sig_smooth = interp1(fas_sig.freq(index_sig_min:index_sig),fas_sig.smooth(index_sig_min:index_sig),xq);
     fas_noi_smooth = interp1(fas_noi.freq(index_noi_min:index_noi),fas_noi.smooth(index_noi_min:index_noi),xq);
     [~, idmin] =  min(abs(fas_sig_smooth - fas_noi_smooth));

     Php = xq(idmin);  

    if ~isempty(Php)
        hp_freq = Php(1,1);
        fprintf('cornerFreqs: high-pass frequency (Hz).... %2.2f\n',hp_freq);
    else     % if no intersection found
        fprintf('cornerFreqs: no high-pass corner frequency detected use 0.1 Hz\n');
        hp_freq = 0.1;
    end

end
%% Compute SNR 
aps = trapz(fas_sig.norm,fas_sig.freq);
apn = trapz(fas_noi.norm,fas_noi.freq);
SNR = 10*log10(aps/apn); %
%% Plotting
if strcmp(debug,'True')
    if isempty(Php)
        hp_freq = 0.1;
    end
    if  lp_freq > 20
        lp_freq = 20;
    end
    fig = figure; clf;
    fsz = 12;
    set(gcf,'position',[300 83 600 500]);
    set(gca,'TickLength',[.0025 .0025]);
    axis square;
    h1 = loglog(fas_sig.freq,fas_sig.norm,'r'); hold on;
    h2 = loglog(fas_noi.freq,fas_noi.norm,'m');
    h3 = loglog(fas_sig.freq,fas_sig.smooth,'k','LineWidth',2);
    h4 = loglog(fas_noi.freq,fas_noi.smooth,'g','LineWidth',2);
    grid on;
    set(gca,'fontname','times','fontsize',fsz);
    xlabel('Frequency, Hz','FontSize',(fsz+1),'fontname','times');
    ylabel('Fourier Amplitude Spectrum, cm/s^2.s','FontSize',(fsz+1),...
        'fontname','times');
    h5 = line([lp_freq lp_freq],ylim,'Color','k','LineStyle','--');
    h6 = line([hp_freq hp_freq],ylim,'Color','b','LineStyle','--');
    legend([h1,h2,h3,h4,h5,h6],('Signal'),('Noise'),('Signal smooth'),...
        ('Noise smooth'),('Low-pass corner'),('High-pass corner'),...
        'Location','northwest');
    titl = strcat('High-pass:', num2str(hp_freq,'%2.2f'),...
        ' Hz & Low-pass: ',num2str(lp_freq,'%2.0f') ,' Hz');
    %axis([0.01, 100, 1e-8 0.01]);
    xlim([0.01, 100])
    title(titl);
    savePDF(plot_path,strcat(plot_name,'_FAS.png'));
   % close(fig)
end
end

%% Compute FAS
function [fas] = computeFAS(w,dt)
%     [N,~]=size(w);
%     
%     ffta=abs(fft(w,N));
%     
%     % first half of the Fourier amplitude
%     ffta_h=ffta(1:floor(N/2));
%     
%     % normalized Fourier amplitude spectrum
%     fas.norm = (ffta_h./N);
%     
%     % frequency vector
%     fas.freq =1/(N*dt)*(1:N/2);

    %% Calculation
    % Nyquist frequency (highest frequency)
    Ny = (1/dt)/2; xgtt = w;
    % number of points in xgtt
    L  = length(xgtt); 
    % Next power of 2 from length of xgtt
    NFFT = 2^nextpow2(L);
    % frequency spacing
    df = 1/(NFFT*dt);
    % Fourier amplitudes 
    U = abs(fft(xgtt,NFFT))*dt; 
    % Single sided Fourier amplitude spectrum
    U = U(2:Ny/df+1);
    % frequency range
    fas.freq = linspace(df,Ny,Ny/df); 
    fas.norm = U;


end

%% Get intersection of two curves
% Adapted and slightly modified from curve intersections - InterX.m written
% by NS, available at
% https://www.mathworks.com/matlabcentral/fileexchange/22441-curve-intersections

function P = getInterX(x1,y1,x2,y2)
% transpose arrays
x1  = x1';
y1  = y1';

% differentiate arrays
dx1 = diff(x1);
dy1 = diff(y1);

dx2 = diff(x2);
dy2 = diff(y2);

% determine signed distances
S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);

% obtain segments where an intersection is expected
arg = dx1.*y2-dy1.*x2;
arg1 = arg(:,1:end-1)-S1;
arg2 = arg(:,2:end)-S1;
C1 = feval(@le,arg1.*arg2,0);

clear arg;
S2 = S2';
arg = (y1.*dx2-x1.*dy2)';
arg1 = arg(:,1:end-1)-S2;
arg2 = arg(:,2:end)-S2;
C2 = feval(@le,arg1.*arg2,0);

% find logical AND of two matrices (C1 and C2); the result contains
% logical 1 (true) only where both matrices contain nonzero values
[i,j] = find(C1 & C2'); % i --> rows; j --> cols
if isempty(i), P = zeros(2,0); return; end

% transpose and prepare for output
dx2=dx2';
dy2=dy2';
L = dy2(j).*dx1(i) - dy1(i).*dx2(j);

% solve system of eqs to get the common points, return unique values
P = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
    dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows')';
end

%% Crop and save MatLAB figure as PDF
function savePDF(plot_path,plot_name)
    % check if directory exists, if not create one
    if ~exist(plot_path, 'dir')
        mkdir(plot_path)
    end
    
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(fig,'-dpng',strcat(plot_path,plot_name));
end