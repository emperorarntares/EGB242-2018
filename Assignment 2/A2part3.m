%% Assignment 2, Part 3 (Choosing a landing site)
%  Do not change before line 32
%  You will need to have generated A2P3Data.mat from 
%  GenerateAssignment2Data.m before working with this file.

%  Clearing and preparing the workspace
clear; clc; close all;

%  Load assignment data from A2P3Data.mat.
load('A2P3Data.mat');  

%=================================================================%
%
% Refer to the assignment sheet for details.
% Names of the variables are important,
% e.g. 'a1' is considered a different variable to 'A1'.
% Make sure variables have been declared as they appear in the brief.
%
% t - time domain vector
% k - frequency vector
% SIG - Fourier Transform
% T - selected value
% sigNoise - estimated noise signal
% a0,an,bn - Trigonometric Fourier Series Coefficients
% OR
% c0,cn - Complex Fourier Series Coefficients
% sigNoise_fs - approximation of noise
% im1 - image 1
% im2 - image 2
% 
%====Enter your code below this line================================

%% 2-3.1: View the noisy image

% Display the first image
figure(1)
imshow(reshape(sig(1,:) , [480, 640]));

%% 2-3.2: Reference vectors

% Sampling frequency
fs = 1000; % pixels/sec

% Sampling period
Ts = length(sig)/fs;

% Time vector
t = linspace(0,Ts, length(sig) + 1); t(end) = []; % 307.2 seconds

% Frequency vector
k = linspace(-fs/2,fs/2,length(sig) + 1); k(end) = [];

%% 2-3.3: Visualise the received signal

% Plot the signal
figure(2)

% Time domain
subplot(2,1,1)
plot(t,sig(1,:))
xlim([0 3]) % First 3 seconds
xlabel('Time (s)')
ylabel('Amplitude')
title('Image Signal 1 - Time Domain')

% Frequency domain
SIG = fft(sig(1,:));
subplot(2,1,2)
plot(k,abs(fftshift(SIG))/fs)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Image Signal 1 - Frequency Domain')

%% 2-3.4: Estimate the periodic noise

% Run estimateNoise.p - transpose vector
sigNoise = estimateNoise(sig(1,:),candidateT(1)).';

% Repmat sigNoise approximately around vector length
sigNoise = repmat(sigNoise,[1 220]);

% Remove excess data points
sigNoise(length(sig) + 1:end) = [];

% Period of noise - measured in seconds
T = candidateT(1)/fs;

% Plot to compare sigNoise and received signal
figure(3)

% Signal 1
plot(t,sig(1,:))
xlim([0 3]) % First 3 seconds
hold on

% sigNoise 1
plot(t,sigNoise)
xlim([0 3]) % First 3 seconds
hold off

xlabel('Time (s)')
ylabel('Amplitude')
title('Image Signal and sigNoise 1')
legend('Received signal','sigNoise')

%% 2-3.5: Model the periodic noise

x1 = t(1:candidateT(1)); % Period time vector
noise = sigNoise(1:candidateT(1)); % Signal sample period

f0 = 1/T; % Frequency

% Computing with Trigonometric Fourier series

% Number of harmonics
numHarm = 6;
n_trig = (1:numHarm).';

% Omega value
wnt = 2 * pi * f0 * (n_trig * x1);

% TFS Coefficients
a0 = (1/T) * sum(noise) / fs;
an = (2/T) * noise * cos(wnt).' / fs;
bn = (2/T) * noise * sin(wnt).' / fs;

%% 2-3.6: Bias

% DC offset
a0 = 0;

%% 2-3.7: Generate the approximation

for ii = 1:length(n_trig)
    % Single harmonic component
    xtf = an(ii)*cos(2*pi*f0*n_trig(ii)*x1) + bn(ii)*sin(2*pi*f0*n_trig(ii)*x1);
    % Expanding vector size to t length
    h_component = repmat(xtf,[1 220]);
    h_component(length(sig)+1:end) = [];
    % Inserting harmonics to each row for n = 1 to 6
    n_matrix(ii,:) = h_component;
end

% Approximation of noise
sigNoise_fs = sum(n_matrix);

%% 2-3.8: Compare the approximation

% Plot sigNoise and sigNoise_fs
figure (4)

% sigNoise
plot(t, sigNoise);
xlim([0 3])
hold on

% sigNoise_fs
plot(t, sigNoise_fs);
xlim([0 3])
hold off

xlabel('Time (s)')
ylabel('Amplitude')
title('sigNoise vs sigNoise TFS approximation')
legend('sigNoise', 'sigNoise Approximated')

% The sigNoise has drastically changed when it is approximated 6 times
% shown in this figure. However, the display shows that the altered signal
% has not been approximated enough to fully align against the original
% sigNoise, which requires a high number of harmonics that is more than 6.

%% 2-3.9: De-noise

% Remove noise
im1 = sig(1,:) - sigNoise_fs;

% Display image
figure(5)
imshow(reshape(im1, [480, 640]));

% As a result, the image has seen some changes due to de-noising the image
% with the approximated sigNoise. It is however has not been processed
% significantly to notice a major improvement on the image, but can be
% seen clear enough to tell that it is located around a rocky mountain; more
% thin lines and less darker compared to the original image. To solve this,
% more number of harmonics is needed to calculate a more accurate
% approximation of sigNoise.

% Improved TFS approximation - Number of harmonics: 15
numHarm = 15;
n_trig = (1:numHarm).';

wnt = 2 * pi * f0 * (n_trig * x1);

a0 = 0;
an = (2/T) * noise * cos(wnt).' / fs;
bn = (2/T) * noise * sin(wnt).' / fs;

for ii = 1:length(n_trig)
    xtf = an(ii)*cos(2*pi*f0*n_trig(ii)*x1) + bn(ii)*sin(2*pi*f0*n_trig(ii)*x1);
    h_component = repmat(xtf,[1 220]);
    h_component(length(sig)+1:end) = [];
    n_matrix(ii,:) = h_component;
end

% New approximation of noise
sigNoise_fs = sum(n_matrix);

% Plot sigNoise and sigNoise_fs
figure (6)
plot(t, sigNoise,'r');
xlim([0 3])
hold on
plot(t, sigNoise_fs,'bl');
xlim([0 3])
hold off
xlabel('Time (s)')
ylabel('Amplitude')
title('sigNoise vs sigNoise TFS approximation')
legend('sigNoise', 'sigNoise Approximated')

% Remove noise and display image
im1 = sig(1,:) - sigNoise_fs;
figure(7)
imshow(reshape(im1, [480, 640]));

% The image has seen significant changes now due to increasing the number
% of harmonics to apply for approximation. We can now see the middle of the
% image clearly, where a rocky hill is shown with a mountain background.
% This area is however only shown visibly and both sides of the image appears to
% be grainy, but still seeable.

%% 2-3.10: Remove the bandlimited random noise

% Image 1 FFT
SIG1 = fftshift(fft(im1(1,:)));

% Plot SIG1 to view noise
figure(9)
subplot(2,1,1)
plot(k,abs(SIG1)/fs)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('SIG1 - Bandlimited Random Noise')

% Locating noise manually
Bandnoise_SIG1 = {'Negative Hz', find(k > -256.5 & k < -240.7),...
                  'Positive Hz', find(k < 256.5 & k > 240.7)};

% Cut off noise range
SIG1(74797:79671) = 0; % Negative frequency
SIG1(227531:232405) = 0; % Positive frequency

% SIG1 without noise
subplot(2,1,2)
plot(k,abs(SIG1)/fs)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('SIG1 - Noise Removed')

% Bring SIG1 back to time domain via inverse Fourier Transform
im2(1,:) = ifft(ifftshift(SIG1));

%% 2-3.11: Choose a site

% De-noise the rest of the signal
for ii = 2:length(sig(:,1))
   im1(ii,:) = sig(ii,:) - sigNoise_fs; 
end

% Repeat process from 3.10

% FFT remaining signals
for ii = 2:length(im1(:,1))
    fft_signal = fftshift(fft(im1(ii,:)));
    SIGS(ii,:) = fft_signal;
end

% The matrix to hold 3 Signals
SIGS(1,:) = [];

% Plot each SIG to inspect noise
figure(10)
for ii = 1:length(SIGS(:,1))
    subplot(3,1,ii);
    plot(k,abs(SIGS(ii,:))/fs)
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title(['SIG' num2str(ii) + 1 ' - Bandlimited Random Noise'])
end

% Extract SIGS as individual vectors to modify
SIG2 = SIGS(1,:);
SIG3 = SIGS(2,:);
SIG4 = SIGS(3,:);

% Locate each SIG noise and organise
Bandnoise_SIG2 = {'Negative Hz',find(k > -264.6 & k < -249.2),... 
                  'Positive Hz',find(k < 264.6 & k > 249.2)};
              
Bandnoise_SIG3 = {'Negative Hz',find(k > -246.9 & k < -232),... 
                  'Positive Hz',find(k < 246.9 & k > 232)};

Bandnoise_SIG4 = {'Negative Hz',find(k > -267.2 & k < -251.1),... 
                  'Positive Hz',find(k < 267.2 & k > 251.1)};

Bandnoises = struct('SIG2',(Bandnoise_SIG2),'SIG3',(Bandnoise_SIG3),...
                    'SIG4',(Bandnoise_SIG4));

% Remove noise from frequency
SIG2(72305:77032) = 0; SIG2(230170:234897) = 0;
SIG3(77760:82334) = 0; SIG3(224868:229442) = 0;
SIG4(71518:76448) = 0; SIG4(230754:235684) = 0;

% Allocate SIGS back to their respective rows
SIGS(1,:) = SIG2;
SIGS(2,:) = SIG3;
SIGS(3,:) = SIG4;

% All SIGS without noise
figure(11)
for ii = 1:length(SIGS(:,1))
    subplot(3,1,ii);
    plot(k,abs(SIGS(ii,:))/fs)
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title(['SIG' num2str(ii) + 1 ' - Noise Removed'])
end

% Convert SIGS to time domain
for ii = 1:length(SIGS(:,1))
    ifft_signal = ifft(ifftshift(SIGS(ii,:)));
    im2(ii + 1,:) = ifft_signal;
end

% Display 4 landing sites
figure(12)
for ii = 1:length(im2(:,1))
    subplot(2,2,ii)
    imshow(reshape(im2(ii,:), [480, 640]))
    title({'Landing Site '; ii})
end

% Landing site 2 looks like the best location to land due to smooth
% terrain, whereas landing site 1 is very rocky and landsites 3 & 4 have
% mountains which may prove to be dangerous and unstable for landing.

%% 2-3.12: Resolution

% Plotting the images individually to be able to read navigational numbers
for ii = 1:length(im2(:,ii))
    figure()
    imshow(reshape(im2(ii,:), [480, 640]))
    title(['Landing Site ' num2str(ii)])
end

% Image 1: +946 , +219 , +656 , +184
% Image 2: +445 , +422 , +136 , +979
% Image 3: +605 , +383 , +612 , +516
% Image 4: +246 , +407 , +084 , +163
