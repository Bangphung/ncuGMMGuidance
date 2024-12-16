function [x_f] = bandpass(c, flp, fhi, dt, n)
        %  Butterworth Acausal Bandpass Filter
        %
        %   Syntax:
        %          x_f = bandpass(c, flp, fhi, dt, n)
        %
        %   Input:
        %            x = input time series
        %          flp = low-pass corner frequency in Hz
        %          fhi = high-pass corner frequency in Hz
        %           dt = sampling interval in second
        %            n = order
        %
        %   Output:
        %        [x_f] = bandpass filtered signal
        fnq = 1/(2*dt);              % Nyquist frequency
        Wn = [flp/fnq fhi/fnq];      % Butterworth bandpass non-dimensional frequency
        [b, a] = butter(n,Wn);
        x_f = filtfilt(b,a,c);
return