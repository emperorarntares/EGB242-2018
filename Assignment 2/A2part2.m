%% Assignment 2, Part 2 (EEG signal analysis)
%  Do not change before line 35
%  You will need to have generated A2P2Data.mat from 
%  GenerateAssignment2Data.m before working with this file.

%  Clearing and preparing the workspace
clear; clc; close all;

%  Load assignment data from A2P1Data.mat.
load('A2P2Data.mat');  

%=================================================================%
%
% Refer to the assignment sheet for details.
% Names of the variables are important,
% e.g. 'a1' is considered a different variable to 'A1'.
% Make sure variables have been declared as they appear in the brief.
%
% Ts - sampling period
% t - time domain vector
% MUX - Fourier Transform
% k - frequency vector
% fshift - frequency shifts
% Mag - magnitude shifts
% Phishift - phase shifts
% xdm - EEG data
% XDM - Fourier Transform
% freqResponse - frequency response of systems
% imp - impulse response of chosen system
% EEG - filtered signals
% eeg - time domain of filtered signals
% eqConv - convolution
%
%====Enter your code below this line================================

%% 2-2.1: Represent the spectrum analyzer

% Insert muxSignal as row vector
muxSignal = muxSignal';

% muxSignal points
samples = length(muxSignal);

% Sampling period
Ts = 1/fs;

% Time domain vector
t = linspace(0,samples*Ts,samples+1); t(end) = [];

% Plot muxSignal
figure(1)
plot(t, muxSignal);
xlabel('Time(s)');
ylabel('Amplitude');
title('muxSignal');

% Compute FFT of muxSignal
MUX = fft(muxSignal);
MUX = fftshift(MUX)/fs;

% Frequency vector
k = linspace(-fs/2,fs/2,samples+1); k(end) = [];

% Plot MUX
figure(2)
plot(k, abs(MUX));
xlabel('Frequency(Hz)');
ylabel('Amplitude');
title('Magnitude spectrum of MUX');

%% 2-2.2: Determine demultiplexing parameters

% Find peaks and frequency shifts of MUX
[peaks, location] = findpeaks(abs(MUX),k, 'NPeaks', 10, 'MinPeakHeight', 0.08);

% Frequency shifts
fshift = location(location>=0); % Positive values
fshift = sort(fshift); % Ascending order

% Magnitude shifts
Mag = peaks(6:end);

% Phase shifts
Phishift = angle(MUX(location>=0));

%% 2-2.3: Remove the frequency shifts for all five astronauts

% Call FDMDemux.p
xdm = FDMDemux(muxSignal,t,Mag,fshift,Phishift);

%% 2-2.4: Review

% Compute FFT in xdm
for ii = 1:length(fshift)
    XDM(ii,:) =  fftshift(fft(xdm(ii,:)));
end

% Plot Magnitude Spectrum for each Data Stream
figure(3)
for i = 1:length(fshift)
    subplot(3,2,i);
    plot(k,abs(XDM(i,:))/fs);
    xlabel('Frequency(Hz)');
    ylabel('Magnitude');
    title(['Data Stream Magnitude ' num2str(i)]);
end

%% 2-2.5: Mathematical analysis

% Display transfer function values
for i = 1:length(sys)
   disp(i);
   factorTF(sys(i));
end

% Comment on each system

%% 2-2.6: System Analysis

% Impulse responses
figure(4)
for i = 1:length(sys)
   subplot(2,2,i);
   impulse(sys(i));
   title(['Impulse response of System ' num2str(i)]);
end

% Step responses
figure(5)
for i = 1:length(sys)
   subplot(2,2,i);
   step(sys(i));
   title(['Step response of System ' num2str(i)]);
end

% Bode plots
figure(6)
for i = 1:length(sys)
   subplot(2,2,i);
   bode(sys(i));
   title(['Bode diagram of System ' num2str(i)]);
end

% Pole zero map
figure(7)
for i = 1:length(sys)
   subplot(2,2,i);
   pzmap(sys(i));
   title(['Pole-zero map of System ' num2str(i)]);
end

%% 2-2.7: Recommend

% See report

%% 2-2.8: Filter the signals

% Setup frequency response for our chosen system

% Complex variable
s = 1i*2*pi*k; % k = frequency vector

% System 3
freqResponse = (413829.0985)./((s + 643.2955).*(s + 643.2955));

% Multiply XDM with freqResponse and store in EEG matrix
for ii = 1:length(XDM(:,1))
   EEG(ii,:) = freqResponse .* XDM(ii,:); 
end

% Convert EEG back to time domain
for ii = 1:length(EEG(:,1))
   eeg(ii,:) = real(ifft(ifftshift(EEG(ii,:)).*fs)); 
end

% Plot each EEG alongside with its eeg
for ii = 1:length(EEG(:,1))
    figure()
    
    % EEG frequency domain
    subplot(2,1,1)
    plot(k, abs(EEG(ii,:))/fs) 
    xlabel('Frequency(Hz)')
    ylabel('Magnitude')
    title(['EEG Signal ' num2str(ii) ' (Frequency Domain)']);
    xlim([-50, 50])
    
    % eeg time domain
    subplot(2,1,2)
    plot(t, eeg(ii,:)/fs) 
    xlabel('Time(s)')
    ylabel('Amplitude')
    title(['eeg Signal ' num2str(ii) ' (Time Domain)']);
end

%% 2-2.9: Equivalence with convolution

% Impulse response
imp = impulse(sys(3),t);

% Convolution of astronaut 1
eqConv = conv(imp,xdm(1,:))/fs;
eqConv = eqConv(1:samples); % Same length

% Compare eqConv with eeg
figure(13)

% eqConv
subplot(2,1,1)
plot(t,eqConv)
xlabel('Time(s)')
ylabel('Amplitude')
title('eqConv of Astronaut 1 (Time Domain)')

% eeg
subplot(2,1,2)
plot(t,eeg(1,:)/fs)
xlabel('Time(s)')
ylabel('Amplitude')
title('eeg Signal of Astronaut 1')

% The plot shows no difference between using convolution and multiplying
% the frequency with time domain. It states that multiplication in the time
% domain also corresponds to convolution in the frequency domain, hence it
% proves so to be true.

%% 2-2.10: Compare

% See report

%% 2-2.11: Digital de-noising

% Plotting the whole EEG in one figure
figure(14)
for ii = 1:length(EEG(:,ii))
    plot(k, abs(EEG(ii,:))/fs) 
    hold on
end
hold off
xlim([-50 50]) % Frequency range
xlabel('Frequency(Hz)') 
ylabel('Magnitude')
title('EEG Astronaut Signals') 
for i = 1:5
    astroNum{i} = ['Astronaut ' num2str(i)];
end
legend(astroNum)

% Beyond bandlimit of 35Hz elements
noise = find(k < -35 | k > 35);

% Creating a filter for corrupted signal
SIG4 = EEG(4,:);

% Removing noise identified outside of bandwidth
SIG4(noise) = 0;

% Restore signal back to EEG
EEG(4,:) = SIG4;

% Plot EEG4
figure(15)
plot(k,abs(EEG(4,:))/fs)
xlim([-50 50]) % Frequency range
xlabel('Frequency(Hz)') 
ylabel('Magnitude')
title('Signal 4 (filtered)') 
legend('Astronaut 4')

% Filtered eeg of Signal 4
eeg(4,:) = real(ifft(ifftshift(EEG(4,:).*fs))); 

% Final representation of Astronaut 4
figure(16)

% EEG frequency domain
subplot(2,1,1)
plot(k, abs(EEG(4,:))/fs) 
xlabel('Frequency(Hz)')
ylabel('Magnitude')
title('EEG Signal 4 (Frequency Domain)');
xlim([-50, 50])
    
% eeg time domain
subplot(2,1,2)
plot(t, eeg(4,:)/fs) 
xlabel('Time(s)')
ylabel('Amplitude')
title('eeg Signal 4 (Time Domain)');


%% 2-2.12: Visual analysis

% See report
