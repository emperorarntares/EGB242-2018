%% Assignment 1 Part A - Section A3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%
%  Do not change before line 28.
%  If you have not generated Data1A from GenerateDataAssignment1A.m,
%  do that now.

%  Clearing and preparing the workspace
clear; clc; close all;

%  Load assignment data.
load Data1A;  

% VARIABLES:
% t - Time vector
% T - Period
% additive noise - Your noise waveform
% a0, an, bn - Trig Fourier series variables
% OR
% c0, cn - Complex Fourier series variables
% FS1 - Fourier series approximation
% dnSnd - De-noised resulting wave

%==================================================================
% Refer to the assignment sheet for details on variable naming.
% Names of the variables are important,
% e.g. 'a1' is considered a different variable to 'A1'.
%====Enter your code below this line================================

%% A3.0 - Setting up common variables for mission.m
T = 5; %Period
f0 = 1/T; %frequency
samples = 44100; %rate per second

x = linspace(0, T, 5*samples + 1); x(end) = []; %Time vector (5s)
t = linspace(0, 4*T, 20*samples + 1); t(end) = []; % t - Time vector (20s)

Ts = t(2) - t(1); %sampling frequency

%% A3.1 Identify noise waveform using a test signal

%plot noiseSound
figure;
plot(t, noiseSound)
hold on

%set up s3
s3(x<=2.5) = C;
s3(x>2.5) = 6;
s3_hinf = repmat(s3,[1 4]); %cycles 4 times

%plotting s3 against noiseSound
plot(t,s3_hinf,'--','LineWidth',2)
xlabel('Time(s)')
ylabel('Amplitude')
title('noiseSound vs S_3(t) Waveform');
legend('noiseSound','Test Signal (S3)');

%% A3.2 Generate noise waveform
additive_noise = noiseSound;

%% A3.3 - Using Complex Fourier series to evaluate noise signal
numHarm = 10;
%setting up trigonometric coefficients for ESF
n_trig = (1:numHarm).';
a0 = (1/T)*sum(s3)*Ts;
an = (2/T)*s3*cos(2*pi*f0*n_trig*x).'*Ts;
bn = (2/T)*s3*sin(2*pi*f0*n_trig*x).'*Ts;

%apply variables to complex coefficients
c0 = a0;
cn_pos = 1/2*(an - 1i*bn);
cn_neg = 1/2*(an + 1i*bn);
n_comp = (-numHarm:numHarm).';
cn = [fliplr(cn_neg), c0, cn_pos];

%% A3.4 - Fourier series approximation (FS1)
for ii = 1:length(n_comp)
    fs1_app(ii,:) = cn(ii) * exp(1j * 2 * pi * f0 * n_comp(ii) * t);
end

FS1 = sum(fs1_app); %signal noise approximation

%plotting Fourier series approximation
figure;
plot(t,additive_noise)
hold on
plot(t,real(FS1))
xlabel('Time(s)')
ylabel('Amplitude')
title('Fourier Series Approximation')
legend('noiseSound','FS1');

%% A3.5 - Recover corrupted speech
%Transpose to minimise matlab resources
additive_noise = additive_noise'; 
%de-noised result
dnSnd = additive_noise-FS1;

%% A3.6 - Plot and listen to recovered Speech Signal
figure;
plot(t,real(dnSnd),'g')
xlabel('Duration(s)')
ylabel('Speech')
title('De-noised Signal')

%play audio
soundsc(real(double(dnSnd)), 44100)

