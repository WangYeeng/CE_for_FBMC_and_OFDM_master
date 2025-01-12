%==========================================================================
% Author: Ying Wang
% Date: 2023.12.18 
% Email: wangyingstu@163.com
% Function:For checking the orthogonality of different systems.
%==========================================================================
clear; clc; close all;
addpath('./MultiCarrierSystem',...
        './ChannelEstimation','./Channel');
%% Parameters
L                   = 25;                           % Number of subcarriers
SubcarrierSpacing   = 15e3;                         % Subcarrier spacing (15kHz, same as LTE)
CarrierFrequency    = 2.5e9;                        % Carrier frequency
K                   = 8;                           % Number of symbols for FBMC
K_OFDM              = 6;                        % Number of symbols for OFDM (same as in LTE)
CP_Length           = 1/SubcarrierSpacing/14;       % LTE CP Length (seconds)
CP_Length_FBMC_OTFS = 2;                            % The CP length for FBMC-based OTFS must be a power of two.
SamplingRate        = 15e3*14*12*5;                 % Sample rate 

%% Channel (different channel environments)
Velocity_kmh        = 500;                          % Velocity (km/h)
PowerDelayProfile   = 'ExtendedVehicularA';         % Power Delay Profile
NrWSSUS             = 200;

%% FBMC 
ConFBMC = FBMC(...
    L,...                                           % Number of subcarriers
    K,...                                           % Number of symbols for FBMC
    SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
    SamplingRate,...                                % Sample rate (Samples/s)
    0,...                                           % Intermediate frequency of the first subcarrier (Hz)
    false,...                                       % Transmit real valued signal
    'PHYDYAS-OQAM',...                              % Prototype filter
    4, ...                                          % Overlapping factor
    0, ...                                          % Initial phase shift
    true ...                                        % Efficient IFFT implementation, true or false
    );

ZeroGuardTimeLength = ((ConFBMC.Nr.SamplesTotal-(round(SamplingRate/SubcarrierSpacing)+CP_Length*SamplingRate)*K_OFDM)/2)/SamplingRate;
OFDM = OFDM(...
    L,...                                           % Number of available subcarriers
    K_OFDM,...                                      % Number of symbols for OFDM
    SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
    SamplingRate,...                                % Sample rate (Samples/s)
    0,...                                           % Intermediate frequency of the first subcarrier (Hz)
    false,...                                       % Transmit real valued signal
    CP_Length, ...                                  % Cycle prefix length (s) 1/SubcarrierSpacing/(K/2-1)
    ZeroGuardTimeLength ...                         % Zero guard length (s)
    );

N = OFDM.Nr.SamplesTotal;

%% Channel
ChannelModel = FastFading(...
         SamplingRate,...                                   % Sample rate (Samples/s)
         PowerDelayProfile,...                              % Power offset profile
         N,...                                              % Total sample number
         Velocity_kmh/3.6*CarrierFrequency/2.998e8,...      % Maximum Doppler shift: Velocity_kmh/3.6*CarrierFrequency/2.998e8
         'Jakes',...                                        % Doppler model: 'Jakes', 'Uniform', 'Discrete-Jakes', 'Discrete-Uniform'. For "Discrete-", we assume a discrete Doppler spectrum to improve the simulation time. This only works accuratly if the number of samples and the velocity is sufficiently large
         NrWSSUS,...                                        % Number of paths for the WSSUS process
         1,...                                              % Number of transmit antennas
         1,...                                              % Number of receive antennas
         0 ...                                              % Gives a warning if the predefined delay taps of the channel do not fit the sampling rate. This is usually not much of a problem if they are approximatly the same.
         );
H = ChannelModel.GetConvolutionMatrix{1};

%% Get transmitted and received matrices
GRX_OFDM = OFDM.GetRXMatrix';
GTX_OFDM = OFDM.GetTXMatrix;
G_FBMC = ConFBMC.GetTXMatrix;

%% Check orthogonality
H_OFDM = GRX_OFDM'*H*GTX_OFDM;
G_FBMC = G_FBMC/(sqrt(1*G_FBMC(:,1)'*G_FBMC(:,1)));
H_FBMC = G_FBMC'*H*G_FBMC;

%% Plotting
Threshold = -60;
figure(1);
subplot(1, 2, 1);
imagesc(10*log10(abs(H_OFDM).^2),[Threshold 0]);
c = colorbar('Ticks',[-60,-50,-40,-30,-20,-10,0]);
c.Label.String = 'dB';
title('OFDM');
xlabel('column index $\ell $','Interpreter','latex');
ylabel('row index $\ell $','Interpreter','latex');
set(gca,"YTick",[1 50 100 150 200 250 300],'XTick',[1 50 100 150 200 250 300],'FontName','Times New Roman','FontSize',13);
subplot(1, 2, 2);
imagesc(10*log10(abs(H_FBMC).^2),[Threshold 0]);
c = colorbar('Ticks',[-60,-50,-40,-30,-20,-10,0]);
c.Label.String = 'dB';
title('FBMC');
xlabel('column index $\ell $','Interpreter','latex');
ylabel('row index $\ell $','Interpreter','latex');
set(gca,"YTick",[1 50 100 150 200 250 300],'XTick',[1 50 100 150 200 250 300],'FontName','Times New Roman','FontSize',13);
