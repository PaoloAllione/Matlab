%%%%%%%%%%%% Applied Signal Processing Laboratory %%%%%%%%%%%%%%%%%%%%%%%%%

% Written by Oliviero Vouch
% Dept. of Electronics and Telecommunications
% Politecnico di Torino
% 2024

% Exercise 1 - T/F Analysis and processing of sinusoidal signals (truncated)


%% Clean workspace
clearvars        % Deletes all the variables in the current workspace
close all force  % Closes all the MATLAB windows except for the IDE     
clc              % Clean the Command Window (but not the hystory)

%% Define simulation parameters
Tmin = 0; % Lower Time limit [sec]
Tmax = 3; % Upper Time limit [sec]
fs=1e4;  % sampling frequency - [samples/second]
Ts=1/fs;  % sampling time (i.e., time resolution) - [seconds]
N=Tmax*fs; % Block size (total number of samples)


%% Task 1: create time axis (according to simulation parameters) 
t=Tmin:Ts:Tmax-Ts;   %time axis

%% Task 2: design a (truncated) sinusoidal signal on the defined time axis

A1 = 2;    % Defines Amplitude
f1 = 200;  % Define the fundamental frequency [Hz] -> how many signal periods I have in each second
Np1 = Tmax*f1; % Integer number of periods inside the window of signal observation

% --- Sinusoidal signal
s1 = A1 * cos( 2 * pi * f1 * t );

%% Task 3: plot the signal in the time-domain and add details
figure('Name','Signal plot [Time domain]') % creates an empty figure with title
plot(t,s1,'b','LineWidth',2)               % creates a Cartesian plot using time axis as x-axis and signal values as y-axis 

% Retrieve hanlde of Cartesian axes in current figure
ax = gca; %gap current axis

% Set labels to the current axes
xlabel(ax,'Time [s]')
ylabel(ax,'s_1(t)')

% Set title to the current axes
title(ax,sprintf('s_1 = A_1cos(2\\pif_1t) ; A_1 = %d, f_1 = %d',A1,f1))

% Display axes grid lines 
grid(ax,"minor")

% Highlight x- and y-axis
xL = xlim;
yL = ylim;
line([0 0], yL,'Color','black');  %x-axis
line(xL, [0 0],'Color','black');  %y-axis

% Set axes fontsize
ax.FontSize = 16;

%% Subtask 3.1: visualize 5 periods of the signal 
ax.XLim = [Tmin 5/f1];

%% Subtask 3.2: use "stem" to visualize discrete-time signal
hold on
stem(t,s1,"filled",'Color','r')

%% Task 4: compute the energy of the signal (compare to nominal)

E_true = A1^2*Tmax/2; % energy of truncated sinusoidal signal

E_mdp = Ts*sum(s1.^2); % midpoint-rule on the time-axis
E_tpz = (Ts/2)*(s1(1)^2+s1(end)^2+2*sum(s1(2:end-1).^2)); % trapezoidal rule

%% Task 5: generate and visualize a C-major chord in octave 5 : composto da tre toni, ovvero tre segnali sinusoidali

% --- Set amplitudes of the three fundamental tones (major triad)
Ac = 1;
Ae = 0.5;
Ag = 0.3;

% --- Notes frequencies can be retrieved from the web
% --- Example of accurate website: https://pages.mtu.edu/~suits/notefreqs.html
f_C5 = 523.25;
f_E5 = 659.25;
f_G5 = 783.99;
% da considerare come frequenze date

%% Subtask 5.1: generate the fundamental tones: creiamo una matrice 
% per ogni riga la matrice contiene i valori di un segnale 
s_tone(1,:) = Ac * cos( 2 * pi * f_C5 * t);
s_tone(2,:) = Ae * cos( 2 * pi * f_E5 * t);
s_tone(3,:) = Ag * cos( 2 * pi * f_G5 * t);

%% Subtask 5.2: generate the chord

% Using for loop (not recommended)
s_Cmaj = zeros(1,length(t));
for i = 1:height(s_tone)
    s_Cmaj = s_Cmaj + s_tone(i,:);
end
% sum matrix columns
s_Cmaj = sum(s_tone,1); % stiamo sommando le colonne della matrice.
% questo ci permette si sommare per ogni istante i valore dei tre segnali 

%% Subtask 5.3: compute the energy of the chord
E_Cmaj_time = (Ts/2)*(s_Cmaj(1)^2+s_Cmaj(end)^2+2*sum(s_Cmaj(2:end-1).^2)); % trapezoidal rule

%% Subtask 5.4: Plot the chord using discrete-time sequence visualization
figure('Name','C5-major chord [Time domain]') % creates an empty figure with title
stem(t,s_Cmaj,'filled','Color','b')              % creates a Cartesian plot using time axis as x-axis and signal values as y-axis 

% Retrieve hanlde of Cartesian axes in current figure
ax = gca;

% Set labels to the current axes
xlabel(ax,'time [nTs]')
ylabel(ax,'s_{C5-maj}[n]')

% Set title to the current axes
title(ax,sprintf('s_{C5-maj} = A_ccos(2\\pif_{C5}t) + A_ecos(2\\pif_{E5}t) + A_gcos(2\\pif_{G5}t)'))

% Display axes grid lines 
grid(ax,"minor")

% Highlight x- and y-axis
xL = xlim;
yL = ylim;
line([0 0], yL,'Color','black');  %x-axis
line(xL, [0 0],'Color','black');  %y-axis

% Set axes fontsize
ax.FontSize = 16;

%% Subtask 5.5: visualize 2 periods of the signal (!!!AGGIUSTARE!!!)
ax.XLim = [Tmin 20e-3];

%% Task 6: play the "C5-major chord" signal
soundsc(s_Cmaj, fs)

%% Subtask 6.1: is the same as playing them separately?
C5 = s_tone(1,:);
E5 = s_tone(2,:); 
G5 = s_tone(3,:); 

soundsc(C5, fs)
soundsc(E5, fs)
soundsc(G5, fs)

%% Subtask 6.2: plot the three notes (hint: use subplots)

figure('Name','Signal plot | C5-E5-G5 [Time domain]')

% C5 tone
subplot(3,1,1)
stem(t,C5,'filled','Color','b')            
sax1 = gca;
xlabel(sax1,'time [nTs]')
ylabel(sax1,'s_{C5}[n]')
% --- Configure the plot
xL = xlim;
yL = ylim;
line(sax1,[0 0], yL,'Color','black');
line(sax1, xL, [0 0],'Color','black');
grid(sax1,"minor")
sax1.FontSize = 16;
sax1.XLim = [Tmin 20e-3]; % first 20 ms of the signal
title(sax1,sprintf('s_{C5} = A_ccos(2\\pif_{C5}t)'))

% E5 tone
subplot(3,1,2)
stem(t,E5,'filled','Color','b')            
sax2 = gca;
xlabel(sax2,'time [nTs]')
ylabel(sax2,'s_{E5}[n]')
% --- Configure the plot
xL = xlim;
yL = ylim;
line(sax2,[0 0], yL,'Color','black');
line(sax2, xL, [0 0],'Color','black');
grid(sax2,"minor")
sax2.FontSize = 16;
sax2.XLim = [Tmin 20e-3]; % first 20 ms of the signal
title(sax2,sprintf('s_{E5} = A_ecos(2\\pif_{E5}t)'))

subplot(3,1,3)
stem(t,G5,'filled','Color','b')            
sax3 = gca;
xlabel(sax3,'time [nTs]')
ylabel(sax3,'s_{G5}[n]')
% --- Configure the plot
xL = xlim;
yL = ylim;
line(sax3,[0 0], yL,'Color','black'); 
line(sax3, xL, [0 0],'Color','black');
grid(sax3,"minor")
sax3.FontSize = 16;
sax3.XLim = [Tmin 20e-3]; % first 20 ms of the signal
title(sax3,sprintf('s_{G5} = A_gcos(2\\pif_{G5}t)'))

%% Task 7: compute the spectrum of the "C5-major chord" signal (frequency domain analysis)
fres=fs/N;          % frequency spacing/resolution (depends only on the observation time T_{max})
f=0:fres:fs-fres;   % frequency axis
ff=f-fs/2;          % symmetric frequency axis (fundamental interval)

S = fft(s_Cmaj);    % FFT spectrum of the "C5-major chord" signal
M = abs(S*Ts);      % magnitude of the fft spectrum
MM = fftshift(M);   % two-sided amplitude spectrum (even function of frequency) --> show the amplitude spectrum in [-fs/2:fs/2-fres]

%% Subtask 7.1: compute the energy and verify Parseval identity
E_Cmaj_freq = (fres/2)*(M(1)^2+M(end)^2+2*sum(M(2:end-1).^2)); % trapezoidal rule

%% Subtask 7.2: Real signal has even amplitude spectrum --> plot the single-sided amplitude spectrum
M_ss = M(1:N/2+1);
M_ss(2:end-1) = 2*M_ss(2:end-1);
f_ss = fres*(0:1:N/2);

figure('Name','Single-sided amplitude spectrum of C5-major chord [frequency domain]') % creates an empty figure with title
plot(f_ss,M_ss,'r','LineWidth',2)
ax = gca;
xlabel(ax,'frequency [Hz]');
ylabel(ax,'|S_{C5-maj}(f)|');
grid(ax,"minor")
ax.FontSize = 16;
xlim([0 +fs/2])
title(ax,sprintf('Single-sided amplitude spectrum of the "C5-major chord" signal'))

% %% Subtask 7.3: Analyze amplitudes of spectral lines
% E_C5_time = (Ts/2)*(C5(1)^2+C5(end)^2+2*sum(C5(2:end-1).^2)); % trapezoidal rule
% E_E5_time = (Ts/2)*(E5(1)^2+E5(end)^2+2*sum(E5(2:end-1).^2)); % trapezoidal rule
% E_G5_time = (Ts/2)*(G5(1)^2+G5(end)^2+2*sum(G5(2:end-1).^2)); % trapezoidal rule