%% Assignment 1 Part A - Section A2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%
%  Do not change before line 31.
%  If you have not generated Data1A from GenerateDataAssignment1A.m,
%  do that now.

%  Clearing and preparing the workspace
clear; clc; close all;

%  Load assignment data.
load Data1A;  

% VARIABLES:
% t - Time vector
% T - Period
% n_trig - Number of harmonics for Trigonometric Fourier Series
% a0, an, bn - Trig Fourier series variables
% s2_hinf - s2(t) ideal time series representation 
% s2_matrix - s2(t) harmonic component matrix
% s2_approx - s2(t) signal approximation
% n_comp - Number of harmonics for Complex Fourier Series 
% c0, cn - Complex Fourier series variables
% s3_hinf - s3(t) ideal time series representation 
% s3_matrix - s3(t) harmonic component matrix
% s2_approx - s3(t) signal approximation

%==================================================================
% Refer to the assignment sheet for details on variable naming.
% Names of the variables are important,
% e.g. 'a1' is considered a different variable to 'A1'.
%====Enter your code below this line================================

%% A2.0 - Setting up common variables for preparation.m
T = 5; %period (t)
f0 = 1/T; %frequency
samples = 1e2; %100 points per period

x = linspace(0,T,samples + 1); x(end) = []; %1 period time vector
t = linspace(0,5*T,5*samples + 1); t(end) = []; %time vector

numHarm = 3; %no. of harmonics for Fourier series approximation
fs = 1/(x(2) - x(1)); %sampling frequency

if numHarm <= 0
   error('Number of Harmonics must be greater than 0.');
end

%% A2.1 - Generate and plot a periodic signal based on s2(t) named s2 hinf.
s2 = exp(-(x-B)/4); %s2(t)

%s2 sample plot
figure;
plot(x,s2,'b')
xlabel('t (seconds)')
ylabel('S_2(t)')
title('S_2(t) Single Period');

%signal period spanning 5 cycles
s2_hinf = repmat(s2,[1 5]);

%periodic signal plot of s2
figure;
plot(t,s2_hinf,'b')
xlabel('t (seconds)')
ylabel('S_2(t)')
title('S_2(t) Waveform');

%% A2.2 - Compute the trigonometric coeffcients of s2(t), numerically using MATLAB.
n_trig = (1:numHarm).'; %no. of harmonics for TFS
wnt = 2 * pi * f0 * (n_trig * x);

%trigonometric coefficients for s2(t)
a0 = (1/T) * sum(s2) / fs;
an = (2/T) * s2 * cos(wnt).' / fs;
bn = (2/T) * s2 * sin(wnt).' / fs;

%% A2.3 - Create a 4x500 matrix called s2_matrix.
for ii = 1:length(n_trig)
    xtf = an(ii)*cos(2*pi*f0*n_trig(ii)*x) + bn(ii)*sin(2*pi*f0*n_trig(ii)*x); %single harmonic component
    h_component = repmat(xtf,[1 5]); %expanding vector to 500 samples
    s2_matrix(ii,:) = h_component; %inserting harmonics to each row for n = 1,2,3
end

row1(1:1,1:500) = a0; %DC component for first row (n = 0)
s2_matrix = [row1; s2_matrix]; %final 4x500 matrix

%% A2.4 - Now create a vector s2 approx which contains an approximation of the signal s2(t).
s2_approx = sum(s2_matrix);

%% A2.5 - Plot the following Fourier series approximation.
figure;
plot(t,s2_hinf, 'b')
hold on

%computing signal approximations
n = a0; %DC component
for ii = 2:length(n_trig)+1 %iterating from 2nd row of matrix -> last row
    n = n + s2_matrix(ii,:); %signal approximation
    s2_n(ii,:) = n; %no. of approximations of s2(t)
end
s2_n(1,:)=[]; %remove empty row
    
%plotting 3 signal approximations
plot(t,s2_n);
legend('S2 Signal','a0+n1','a0+n1+n2','a0+n1+n2+n3')
xlabel('t (seconds)')
ylabel('S_2(t)')
title('S_2(t) TFS Signal Approximations')

%% A2.6 - For s3(t), repeat the steps of A2.1 to A2.5
%s3(t): consider 2 y1 points
s3 = zeros(size(x));
s3(x>= 0 & x<2.5) = -4.5;
s3(x>=2.5 & x<5) = 6;

%s3 sample plot
figure;
plot(x,s3,'b')
axis([0 5 -6 8])
xlabel('t (seconds)')
ylabel('S_3(t)')
title('S_3(t) Single Period');

%signal period spanning 5 cycles
s3_hinf = repmat(s3,[1 5]);

%periodic signal plot of s3
figure;
plot(t,s3_hinf,'b')
axis([0 25 -6 8])
xlabel('t (seconds)')
ylabel('S_3(t)')
title('S_3(t) Waveform');

%trigonometric coefficients for s3(t)
a0 = (1/T) * sum(s3) / fs;
an = (2/T) * s3 * cos(wnt).' / fs;
bn = (2/T) * s3 * sin(wnt).' / fs;

%exponential coefficients for s3(t)
c0 = a0;
cn_pos = 1/2 * (an - 1i * bn);
cn_neg = 1/2 * (an + 1i * bn);

n_comp = (-numHarm:numHarm).';
%[c-3 c-2 c-1 c0 c1 c2 c3];

cn = [fliplr(cn_neg), c0, cn_pos]; %range of cn values

table(n_comp,cn.','VariableNames',{'n','Cn'}); %table view

%creating a 7x500 matrix called s3_matrix
for ii = 1:length(n_comp)
    %complex fourier series
    s3_matrix(ii,:) = cn(ii) * exp(1j * 2 * pi * f0 * n_comp(ii) * t);
end

%signal approximation values for s3
s3_approx = sum(s3_matrix);

%computing signal approximations
for ii = 1:length(n_comp)/2
    harmonics(ii,:) = real(s3_matrix(ii,:) * 2); %single harmonic component
end

n = c0; %DC component
for ii = ((length(n_comp) -1)/2):-1:1 %iterating no. of harmonics -> 1
    n = n + harmonics(ii,:); %signal approximation
    s3_n(ii,:) = n; %no. of approximations of s3(t)
end
s3_n = flipud(s3_n); %fix legend

%plotting Fourier series approximation
figure;
plot(t,s3_hinf,'b')
axis([0 25 -6 8])
hold on

%plotting 3 harmonic sets
plot(t, s3_n(1,:),'--','LineWidth',3); %first signal approximation
for ii = 2:size(s3_n,1)
    plot(t, s3_n(ii,:)) %rest of the approximations
end
legend('S3 Signal','c0+n1','c0+n1+n2','c0+n1+n2+n3')
xlabel('t (seconds)')
ylabel('S_3(t)')
title('S_3(t) EFS Signal Approximations')

%% A2.7

% The number of harmonics applied increases the approximation accuracy of
% the signal. The current number of harmonics (3) is sufficient enough to
% demonstrate the signals and its approximation, and lowering the number of
% harmonics would not prove so much or display an accurate measurement of
% the signal,creating a broad interpretation of the information as shown in
% s3(t) (Figure 2.6) where the first two signal approximations are the
% same. Raising the number of harmonics increases accuracy and show many
% transformations, however it would further narrow down the waveforms of
% the signal, making it smaller until it turns to the ideal periodic
% series.