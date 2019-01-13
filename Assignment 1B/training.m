%% Assignment 1 - Part B, Section B2 (Coding: Training Exercise)
%
%  M-file: training.m
%  This file is for use in Section B2 (Training Exercise) of Assignment 1B.
%
%  Section B3 (Mars Mission: Combining Signals) will have its own matlab file (mission.m).
%  
%  Preparing your workspace
clear; clc; close all;  
%==================================================================
% Names of variables you will need for marking.
% Refer to the assignment sheet for details.
% Names of the variables are important,
% e.g. 'a1' is considered a different variable to 'A1'.
% Make sure you have declared them as they appear in the list below.
% ---------------------------------------------
% ===== Part 1 Variables =====
%  f1	Frequency vector
%  S1f  Spectrum S_1(f)
%  t1	Time vector
%  s1   Signal s_1(t)
%  Ts1	Sampling period
%  fs1	Sampling frequency
%  k1	Frequency vector
%  S1k	Fourier transform of s1
%  f2   Frequency vector
%  S2f  Fourier transform of s_2(t)
%  t2   Time vector
%  s2	Signal s_2(t)
%  S2k  Fourier transform of s2
%  s3   Signal s_3(t)
%  S3k  Fourier transform of s3
%  S3f  Fourier transform of s_3(t)
% ---------------------------------------------
%====Enter your code below this line================================

%% Setting up common variables
m = 10;
A = 4;
phi = pi;

%% B2.1

% 10000 points
samples = 10000; 

% Frequency vector
f1 = linspace(-50,50,samples+1); f1(end) = [];

% Fourier transform S1(f)
S1f = (exp(-3*(1+1i*2*pi*f1)))./((1+1i*2*pi*f1).^2);

% Plot magnitude and phase of S1(f)
figure;
subplot(2,1,1);
stem(f1,abs(S1f));
xlabel('Frequency');
ylabel('Magnitude');
title('Magnitude and Phase Spectrum of S1(f)');
subplot(2,1,2);
plot(f1,angle(S1f));
xlabel('Frequency');
ylabel('Phase');

%% B2.2

% 4000 points
samples = 4000;

% Time vector
t1 = linspace(0,8,samples+1); t1(end) = [];

% Signal s_1(t)
s1 = (t1-3).*exp(-t1).*heaviside(t1-3);

% Plot s1(t)
figure;
plot(t1,s1)
xlabel('Time(s)');
ylabel('Amplitude')
title('Time domain of S1(t)')

% Sampling period
Ts1 = t1(2) - t1(1);

% Sampling frequency
fs1 = 1/Ts1;

% Frequency vector
k1 = linspace(-fs1/2,fs1/2,length(t1)+1); k1(end) = [];

% Compute discrete Fourier transform
S1k = fft(s1);

%% B2.3

% Plot magnitude and phase of S1k
figure;
subplot(2,1,1);
stem(k1,abs(fftshift(S1k))/fs1);
xlabel('Frequency');
ylabel('Magnitude');
title('Magnitude and Phase Spectrum of S1k(f)');

subplot(2,1,2);
plot(k1,angle(fftshift(S1k)));
xlabel('Frequency');
ylabel('Phase');

%% B2.4

% 20000 points
samples = 20000;

% Frequency vector
f2 = linspace(-50,50,samples+1); f2(end) = [];

% S2(f) component sets
eq1 = (exp(-3*(1+2i*pi*(f2-m))))./(2*(1+2i*pi*(f2-m)).^2); % Equation 1
eq2 = (exp(-3*(1+2i*pi*(f2+m))))./(2*(1+2i*pi*(f2+m)).^2); % Equation 2

% Fourier transform S2(f)
S2f = eq1 + eq2;

% Plot magnitude spectrum
figure;
stem(f2,abs(S2f))
xlabel('Frequency');
ylabel('Magnitude');
title('Magnitude Spectrum of S2(f)');

% 6000 points
samples = 6000;

% Time vector
t2 = linspace(0,8,samples+1); t2(end) = [];

% Signal s_2(t)
s2 = (t2-3).*exp(-t2).*heaviside(t2-3) .* cos(2*pi*m*t2);

% Plot s2(t)
figure;
plot(t2,s2);
xlabel('Time(s)');
ylabel('Amplitude')
title('Time domain of S2(t)')

% Compute discrete Fourier transform
S2k = fft(s2);

%% B2.5

% Sampling period
Ts2 = t2(2) - t2(1);
% Sampling frequency
fs2 = 1/Ts2;
% Frequency vector
k2 = linspace(-fs2/2,fs2/2,length(t2)+1); k2(end) = [];

% Plot S2k vs S2f
figure;
stem(k2,abs(fftshift(S2k))/fs2,'r');
hold on
stem(f2,abs(S2f),'b')
xlabel('Frequency');
ylabel('Magnitude');
title('S2k vs. S2f Magnitude');
legend('S2k','S2f');

%% B2.6

% Signal s_3(t)
s3 = (t1-3).*exp(-t1).*heaviside(t1-3) .* cos(2*pi*m*t1) .* A .* cos(2*pi*m*t1 + phi);

% FFT s3
S3k = fft(s3);

% S3(f) component sets
eq1 = (exp(-3*(1+2i*pi*(f2-2*m))+1i*phi))./(2*(1+2i*pi*(f2-2*m)).^2); 
eq2 = (exp(-3*(1+2i*pi*f2))*(2*cos(phi)))./(2*(1+2i*pi*f2).^2);
eq3 = (exp(-3*(1+2i*pi*(f2+2*m))-1i*phi))./(2*(1+2i*pi*(f2+2*m)).^2); 

% Fourier transform S3(f)
S3f = A/2 * (eq1 + eq2 + eq3);

% Plot S3k & S3f
figure;
stem(k1,abs(fftshift(S3k))/fs1,'r')
hold on
stem(f2,abs(S3f),'b')
xlabel('Frequency');
ylabel('Magnitude');
title('S3k vs. S3f Magnitude');
legend('S3k','S3f');

%% B2.7 - See report