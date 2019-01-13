%% Assignment 1 - Part B, Section B3 (Combining Signals)
%  Do not change before line 40
%  You will need to have generated Data1B.mat from 
%  GenerateDataAssignment1B.m before working with this file.

%  Clearing and preparing the workspace
clear; clc; close all;

%  Load assignment data from Data1B.mat
load('Data1B.mat', 'fs', 'muxSignal');% 

%  The variable loaded are:
%     fs      Sampling frequency for Section B3
%  muxSignal  Multiplexed music signals, for Section B3.

%==================================================================
%
% Names of variables you will need for marking.
% Refer to the assignment sheet for details.
% Names of the variables are important,
% e.g. 'a1' is considered a different variable to 'A1'.
% Make sure variables have been declared as they appear in the list below.
% ---------------------------------------------
% ===== Part 2 =====
%     Ts         Sampling period
%     t          time vector
%     MUX        Fourier transform of muxSignal
% MUX_fftshifted Shifted and scaled Fourier Transform
%     k          Frequency vector
%   fshift       frequency shifts (vector)
%   Ashift       amplitude of peaks (vector)
%   Phishift     phases of peaks (vector)
%     xdm        output of FDMDemux (matrix)
%     XDM        Fourier transform of xdm (matrix)
%     B          Bandwidth
% filteredoutput Filtered Output Signal
% decodedtext    The Decoded Text Output
% ---------------------------------------------
%====Enter your code below this line================================

%% B3.1

% muxSignal points
samples = length(muxSignal);

% Sampling period
Ts = 1/fs;

% Construct a time domain vector
t = linspace(0,samples*Ts,samples+1); t(end) = [];

% Plot muxSignal against t
figure;
plot(t, muxSignal);
xlabel('Time(s)');
ylabel('Amplitude');
title('muxSignal');

%% B3.2

% Fourier transform muxSignal
MUX = fft(muxSignal);

% Frequency vector for MUX
k = linspace(-fs/2,fs/2,length(t)+1); k(end) = [];

% Shifted MUX for plot
MUX_shift = abs(fftshift(MUX)/fs);

% Plot magnitude spectrum
figure;
plot(k,MUX_shift);
xlabel('Frequency');
ylabel('Magnitude');
title('Magnitude Spectrum of MUX');

%% B3.3

% Find peaks of MUX_fftshifted
[Muxpeaks, Muxt] = findpeaks(MUX_shift,k, 'NPeaks', 8, 'MinPeakHeight', 0.4);

% Defining positive frequencies
fshift = Muxt(5:end);
% Defining Positive peaks
Ashift = Muxpeaks(5:end);

%% B3.4

% Finding fshifts in frequency vector k
for i = 1:length(fshift)
    k_Muxt(i) = find(k==fshift(i));
end

% Magnitude values for peaks - should be same as Ashift
Mag = abs(MUX_shift(k_Muxt));

% Phase values for peaks
Phishift = angle(MUX_shift(k_Muxt));

% Plot new figure
figure;
plot(k,MUX_shift);
hold on
plot(fshift,Mag,'ro','Color','g');
xlabel('Frequency');
ylabel('Magnitude');
title('MUX Peaks in fshift Domain');
legend('MUX','Peaks');

%% B3.5 - Frequency shifting module
[xdm] = FDMDemux(muxSignal,t,Mag,fshift,Phishift);

%% B3.6

% Compute FT in xdm
for ii = 1:length(fshift)
    XDM(ii,:) =  fft(xdm(ii,:));
end

% Subplot each XDM
figure;
for i = 1:length(fshift)
   subplot(2,2,i);
   plot(k,abs(XDM(i,:)));
   xlabel('Frequency');
   ylabel('Magnitude');
   title(['XDM Magnitude for Signal' num2str(i)]);
end

%% B3.7

B = 5*10^4; % Bandwidth value

filteredoutput = A1BLPF(xdm,fs,B); % Filter system

%% B3.8

% Signal 1
A1BTextdecode(filteredoutput(1,:),Ts)

% Signal 2
sound(filteredoutput(2,:),fs)

% Signal 3
A1BTextdecode(filteredoutput(3,:),Ts)

% Signal 4
sound(filteredoutput(4,:),fs)

