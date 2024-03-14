% Script to generate simulated signals of a single diffusion weighted image
% (1 excitation) and to process the generated signal in the same way as 
% time-dependent diffusion (TDD) contrast maps are pocessed.

% MiJo 7.12.2023
%% User input
clear all
% switches:
bvalue = 2; % 1 for b1250, 2 for b2500.
roi = 0; % 0 for wm, 1 for tumor 

% example signal levels from actual measurements:
% S0: signal level from b = 0 s/mm^2 acquisition.
% S1: signal with shorter diffusion time.
% S2: signal with longer diffusion time
S0_mean_wm = 627.52; 
S1_b1250_wm = 252.87; S1_b2500_wm = 159.18;
S2_b1250_wm = 254.29; S2_b2500_wm = 160.64;

S0_mean_tumor = 1534.8;
S1_b1250_tumor = 287.38; S1_b2500_tumor = 132.07;
S2_b1250_tumor = 289.39; S2_b2500_tumor = 135.56;

%% Background info
volunteer_S0_WM = 2298.57;
volunteer_noise = 77.59;
factor = volunteer_noise/volunteer_S0_WM;

N = 100000; % number of generated signals
Nd = 6; % number of gradient directions
Nex = 8; % number of excitations per gradient direction

scaled_noise = factor*S0_mean_wm;

if bvalue == 1
    simMean_wm = {S1_b1250_wm, S2_b1250_wm};
    simMean_tumor = {S1_b1250_tumor, S2_b1250_tumor};
elseif bvalue == 2
    simMean_wm = {S1_b2500_wm, S2_b2500_wm};
    simMean_tumor = {S1_b2500_tumor, S2_b2500_tumor};
end

%% Generate noisy signals
temp = zeros(N, Nd, Nex); 
simS = {temp,temp}; % two signals: "short" and "long" diffusion time

rng('shuffle')
for si = 1:2 % two signals
    for d = 1:6 % 6 gr directions
        for ex = 1:8 % 8 excitations per gr dir
            if roi == 0
                simS{si}(:,d,ex) = simMean_wm{si} + scaled_noise.*randn(N,1); % N x of simulated signal
            elseif roi == 1
                simS{si}(:,d,ex) = simMean_tumor{si} + scaled_noise.*randn(N,1);
            end
        end
    end
end

% Generate S0:
if roi == 0
    S0 = S0_mean_wm + scaled_noise.*randn(N,1);
    fprintf('Simulating WHITE MATTER.\n')
elseif roi == 1
    S0 = S0_mean_tumor + scaled_noise.*randn(N,1);
    fprintf('Simulating TUMOR.\n')
end

fprintf('Example SD for S1:\nOriginal: %.4f - ', std(simS{1}(:,1,1)))

%% Signal processing

% Step 1: DW image along 1 gradient is an average over 8 excitations
tempS1 = mean(simS{1}, 3);
tempS2 = mean(simS{2}, 3);
fprintf('Step 1: %.4f - ', std(tempS1(:,1)))

% Step 2: background gradient correction
S1d3 = zeros(N, 3); S2d3 = S1d3;
for d = 1:3
    S1d3(:,d) = sqrt(tempS1(:,d).*tempS1(:,d+3));
    S2d3(:,d) = sqrt(tempS2(:,d).*tempS2(:,d+3));
end
fprintf('Step 2: %.4f - ', std(S1d3(:,1)))

% Step 3: isotropic image (as a geometric mean)
S1 = nthroot(prod(S1d3, 2), 3);
S2 = nthroot(prod(S2d3, 2), 3);
fprintf('Step 3: %.4f - ', std(S1))

% Step 4: b0 normalization
S1n = S1./S0;
S2n = S2./S0;
fprintf('Step 4: %.4f.\n', std(S1n))

% Step 5: subtraction
TDD = S2n - S1n;

% Results
mean_TDD = mean(TDD);
std_TDD = std(TDD);

fprintf('Mean TDD = %.4f, SD TDD = %.4f\n', mean_TDD, std_TDD)
figure, clf
histogram(TDD)
title(sprintf('simulated TDD_{b%d}: distribution of %d signals', bvalue, N))
xlabel('TDD'), ylabel('N')
grid on
xlim([-0.05 0.05])
