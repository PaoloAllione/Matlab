%% Clean workspace
clearvars        % Deletes all the variables in the current workspace
close all force  % Closes all the MATLAB windows except for the IDE     
clc              % Clean the Command Window (but not the hystory)

%% Define simulation parameters for rectangular pulse
A = 2; % rect. pulse amplitude [-]
T = 1; % rect. pulse duration [sec]
Tmin = 0; % Lower limit of observation time-window [sec]
Tmax = 10*T; % Upper limit of observation time-window [sec]
f0 = 1e3; % modulation frequency [Hz]
fs=1e4;  % sampling frequency - [samples/second]
Ts=1/fs;  % sampling time (i.e., time resolution) - [seconds]

%% create Cartesian time/frequency axes (according to simulation parameters) 
t=Tmin:Ts:Tmax-Ts;  % time axis
N = Tmax*fs; % Total number of samples in time axis
NT = T*fs; % Number of samples of rectangular pulse

fres=fs/N;          % frequency spacing/resolution (depends only on the observation time T_{max})
f=0:fres:fs-fres;   % frequency axis
ff=f-fs/2;          % symmetric frequency axis (fundamental interval)

%% Task 2: design the rectangular pulse with amplitude A and duration T on the defined time axis
s = zeros(1,N);
s(1:NT)=A;   

%% Task 3: plot the signal [time-domain]
figure('Name','Signal plot [Time domain]') % creates an empty figure with title
plot(t,s,'b','LineWidth',2)               % creates a Cartesian plot using time axis as x-axis and signal values as y-axis 

% Retrieve hanlde of Cartesian axes in current figure
ax = gca;

% Set labels to the current axes
xlabel(ax,'Time [s]')
ylabel(ax,'s(t)')

% Set title to the current axes
title(ax,sprintf('s = Ap_T(t) ; A = %d, T = %d',A,T))

% Display axes grid lines 
grid(ax,"minor")

% Set axes fontsize
ax.FontSize = 16;

%% Task 4: plot the signal amplitude spectrum [frequency-domain]
S = fft(s);         % FFT spectrum of the rectangular pulse
M = abs(S*Ts);      % magnitude of the fft spectrum
MM = fftshift(M);   % two-sided amplitude spectrum (even function of frequency) --> show the amplitude spectrum in [-fs/2:fs/2-fres]

% % %% Subtask 7.2: Real signal has even amplitude spectrum --> plot the single-sided amplitude spectrum
% % M_ss = M(1:N/2+1);
% % M_ss(2:end-1) = 2*M_ss(2:end-1);
% % f_ss = fres*(0:1:N/2);

figure('Name','Double-sided amplitude spectrum [frequency domain]') 
plot(ff,MM,'r','LineWidth',2)
ax = gca;
xlabel(ax,'frequency [Hz]');
ylabel(ax,'|S(f)|');
grid(ax,"minor")
ax.FontSize = 16;
xlim([-fs/2 +fs/2])
title(ax,sprintf('Double-sided amplitude spectrum'))

%% Task 5: modulate the signal around f0 and plot the modulated spectrum (vs. original spectrum)
smod = s.*exp(1i*2*pi*f0*t);
Smod = fft(smod);         % FFT spectrum of the rectangular pulse
Mmod = abs(Smod*Ts);      % magnitude of the fft spectrum
MMmod = fftshift(Mmod);   % two-sided amplitude spectrum (even function of frequency) --> show the amplitude spectrum in [-fs/2:fs/2-fres]

figure('Name','Double-sided amplitude spectrum [frequency domain] of original vs. modulated signal')

% original signal
subplot(2,1,1)
sax1 = gca;
stem(ff,MM,'filled','Color','r')            
xlabel(sax1,'frequency [Hz]')
ylabel(sax1,'|S(f)|')
% --- Configure the plot
grid(sax1,"minor")
sax1.FontSize = 16;
sax1.XLim = [-fs/2 +fs/2]; 
sax1.YLim = [0 A];
title(sax1,sprintf('amplitude spectrum of the baseband signal'))

% modulated signal
subplot(2,1,2)
sax2 = gca;
stem(ff,MMmod,'filled','Color','r')            
xlabel(sax2,'frequency [Hz]')
ylabel(sax2,'|S_{mod}(f)|')
grid(sax2,"minor")
sax2.FontSize = 16;
sax2.XLim = [-fs/2 +fs/2]; 
sax2.YLim = [0 A];
title(sax2,sprintf('amplitude spectrum of modulated signal'))