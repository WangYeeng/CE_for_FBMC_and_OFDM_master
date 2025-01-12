% =========================================================================   
% (c) 2024 Ying Wang, wangyingstu@163.com
% =========================================================================  
clc; clear; close all;
addpath('./MultiCarrierSystem',...
        './ChannelEstimation','./Channel');
load("R_H_560.mat");
%% Parameters
% Simulation
M_SNR_dB                  = 0:8:40;                 % Signal-to-Noise Ratio in dB
NrRepetitions             = 1500;                   % Number of Monte Carlo repetition (different channel realizations)             
ZeroThresholdSparse       = 10;                     % Set some matrix values, smaller than "10^(-ZeroThresholdSparse)", to zero.

% FBMC and OFDM parameters
L                         = 24;                     % Number of Subcarriers, one resource block consists of 12 subcarriers (and 0.5ms in time)
K                         = 30;
F                         = 15e3;                   % Subcarrier Spacing in Hz
SamplingRate              = F*L*5;                  % Sampling rate in (Samples/s)
CarrierFrequency          = 2.5e9;                  % Carrier Frequency
Cycliclength              = 1/15e3/14;
NrSubframes               = 1;                      % Number of subframes. One subframe requires 1ms for F=15kHz.                             
QAM_ModulationOrder       = 16;                     % QAM signal constellation order, 4, 16, 64, 256, 1024,...

% Channel estimation parameters
NrIterations              = 5;                      % Number of iterations for interference cancellation scheme.
PilotToDataPowerOffset_OFDM    = 2.15;
PilotToDataPowerOffset_FBMC    = 4.3;
% Channel       
PowerDelayProfile         = 'ExtendedVehicularA';   % Channel model: 'VehicularA', 'VehicularB', 'ExtendedVehicularA'

%% OFDM Object (Add zeroes to the OFDM signal so that it fits the FBMC signal)
ZeroGuardTimeLength = 1.15e-4;
OFDM = OFDM(...
    L,...                           % Number Subcarriers
    14*NrSubframes,...              % Number OFDM Symbols
    F,...                           % Subcarrier spacing (Hz)
    SamplingRate,...                % Sampling rate (Samples/s)
    0,...                           % Intermediate frequency first subcarrier (Hz)
    false,...                       % Transmitreal valued signal
    Cycliclength, ...               % Cyclic prefix length (s)
    ZeroGuardTimeLength ...         % Zero guard length (s)
    );

%% FBMC Object
FBMC = FBMC(...
    L,...                               % Number subcarriers
    K*NrSubframes,...                   % Number FBMC symbols
    F,...                               % Subcarrier spacing (Hz)
    SamplingRate,...                    % Sampling rate (Samples/s)
    0,...                               % Intermediate frequency first subcarrier (Hz)
    false,...                           % Transmit real valued signal
    'PHYDYAS-OQAM',...                  % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM, 
    4, ...                              % Overlapping factor (also determines oversampling in the frequency domain)
    0, ...                              % Initial phase shift
    true ...                            % Polyphase implementation
    );

if (OFDM.Nr.SamplesTotal==FBMC.Nr.SamplesTotal)
    N = OFDM.Nr.SamplesTotal;
else
    disp('OFDM and FBMC must have the same sample number!!')
end

%% QAM Object
QAM = SignalConstellation(QAM_ModulationOrder,'QAM');
PAM = SignalConstellation(sqrt(QAM_ModulationOrder),'PAM');
%% Pilot Matrices, 0="Data", 1=Pilot
% OFDM
ChannelEstimation_OFDM = ProfilePilotSymbols(...
    'Diamond',...                           % Pilot pattern
    [...                                    % Matrix that represents the pilot pattern parameters
    OFDM.Nr.Subcarriers,...                 % Number of subcarriers
    6; ...                                  % Pilot spacing in the frequency domain
    OFDM.Nr.MCSymbols,...                   % Number of FBMC/OFDM Symbols
    3.5 ...                                 % Pilot spacing in the time domain
    ],...                                   
    'linear'...                             % Interpolation(Extrapolation) method 'linear','FullAverage'
    );
PilotMatrix_OFDM = ChannelEstimation_OFDM.PilotMatrix;
PilotMatrix_OFDM = repmat(PilotMatrix_OFDM,[1 NrSubframes]);
% FBMC
ChannelEstimation_FBMC = ProfilePilotSymbols(...
    'Diamond',...                           % Pilot pattern
    [...                                    % Matrix that represents the pilot pattern parameters
    FBMC.Nr.Subcarriers,...                 % Number of subcarriers
    6; ...         % Pilot spacing in the frequency domain
    FBMC.Nr.MCSymbols,...                   % Number of FBMC/OFDM Symbols
    8 ...             % Pilot spacing in the time domain
    ],...                                   
    'linear'...                            % Interpolation(Extrapolation) method 'linear','FullAverage'...
    );
PilotMatrix_FBMC =  ChannelEstimation_FBMC.PilotMatrix;

GuardPilotMatrix_FBMC = PilotMatrix_FBMC;
[a,b] = find(PilotMatrix_FBMC);
for i_pilot = 1:length(a)
    GuardPilotMatrix_FBMC(a(i_pilot)+1,b(i_pilot))=-1;
    GuardPilotMatrix_FBMC(a(i_pilot)-1,b(i_pilot))=-1;
    GuardPilotMatrix_FBMC(a(i_pilot),b(i_pilot)+1)=-1;
    GuardPilotMatrix_FBMC(a(i_pilot),b(i_pilot)-1)=-1;
end

%% Cancel Imaginary Interference At Pilot Position Object (Precoding Matrix)
GuardSymbolStructure = ProfileGuardSymbols(...
    'Guard', ...                                        % Cancellation method
    GuardPilotMatrix_FBMC, ...                          % PilotAndAuxiliaryMatrix
    FBMC.GetFBMCMatrix, ...                             % Imaginary interference matrix
    16, ...                                             % Cancel closest interferers
    PilotToDataPowerOffset_FBMC ...                          % Pilot to data power offset
    );

NrPilotSymbols     = sum(PilotMatrix_FBMC(:)==1);
NrDataSymbols_FBMC = sum(PilotMatrix_FBMC(:)==0);
NrGuardSymbols     = GuardSymbolStructure.NrGuardSymbols;
Overhead_FBMC      = (NrPilotSymbols+NrGuardSymbols)/(FBMC.Nr.Subcarriers*FBMC.Nr.MCSymbols);
Rho_FBMC  = GuardSymbolStructure.PilotToDataPowerOffset*GuardSymbolStructure.DataPowerReduction;

%% Cancel Imaginary Interference At Pilot Position Object (Precoding Matrix)
NrPilotSymbols     = sum(PilotMatrix_OFDM(:)==1);
NrDataSymbols_OFDM = sum(PilotMatrix_OFDM(:)==0);

PilotMapping_OFDM = zeros(numel(PilotMatrix_OFDM),numel(PilotMatrix_OFDM));
PilotMapping_OFDM(PilotMatrix_OFDM(:)==1,1:NrPilotSymbols) = sqrt(PilotToDataPowerOffset_OFDM)*eye(NrPilotSymbols);
PilotMapping_OFDM(PilotMatrix_OFDM(:)==0,NrPilotSymbols+1:end) = eye(NrDataSymbols_OFDM);
PilotMapping_OFDM = PilotMapping_OFDM/sqrt(mean(diag(PilotMapping_OFDM*PilotMapping_OFDM')));
DataPowerReduction_OFDM = numel(PilotMatrix_OFDM)/(NrPilotSymbols*PilotToDataPowerOffset_OFDM+NrDataSymbols_OFDM);

Rho_OFDM     = PilotToDataPowerOffset_OFDM*DataPowerReduction_OFDM;
Overhead_OFDM  = (NrPilotSymbols)/(OFDM.Nr.Subcarriers*OFDM.Nr.MCSymbols) + (Cycliclength)/OFDM.PHY.TimeSpacing;
%% Precalculate Transmit and Receive Matrices
G_OFDM = OFDM.GetTXMatrix;
Q_OFDM = (OFDM.GetRXMatrix)';
GP_OFDM = G_OFDM(:,PilotMatrix_OFDM(:)==1);
QP_OFDM = Q_OFDM(:,PilotMatrix_OFDM(:)==1);
G_OFDM_PilotMapping = G_OFDM*PilotMapping_OFDM;

G_FBMC  = FBMC.GetTXMatrix;
Q_FBMC  = (FBMC.GetRXMatrix)';
GP_FBMC = G_FBMC(:,PilotMatrix_FBMC(:)==1);
QP_FBMC = Q_FBMC(:,PilotMatrix_FBMC(:)==1);
G_Gua   = G_FBMC*GuardSymbolStructure.PrecodingMatrix;

disp('Calculate Correlation Matrixe...');
%% Calculate Correlation Matrix of Pilot Estimates for OFDM
% OFDM
[R_hP_est_OFDM,R_hP_est_noInterference_OFDM] = RhPEstimation( ...
    M_SNR_dB, ...        % Signal-to-Noise Ratio (dB)
    F*L, ...             % Bandwidth
    Rho_OFDM, ...      % Pilot-to-data power offset
    GP_OFDM, ...         % Transmit Matrix at the Pilot
    QP_OFDM, ...         % Receive matrix at the pilot
    G_OFDM_PilotMapping, ...% Transmit Matrix
    R_H, ...          % Channel Correlation Matrix
    SamplingRate, ...    % Sampling Rate
    NrPilotSymbols, ...  % Pilot symbol number
    N ...                % Total number of samples
    );
% FBMC
[R_hP_est_FBMC,R_hP_est_noInterference_FBMC] = RhPEstimation( ...
    M_SNR_dB, ...        % Signal-to-Noise Ratio (dB)
    F*L, ...             % Bandwidth
    Rho_FBMC, ...      % Pilot-to-data power offset
    GP_FBMC, ...         % Transmit Matrix at the Pilot
    QP_FBMC, ...         % Receive matrix at the pilot
    G_Gua, ...           % Transmit Matrix
    R_H, ...          % Channel Correlation Matrix
    SamplingRate, ...    % Sampling Rate
    NrPilotSymbols, ...  % Pilot symbol number
    N ...                % Total number of samples
    );

disp('Calculate induction factor Solution and that of no Interference...');
%% Calculate MMSE Estimation Matrix and that of no Interference 
% OFDM
[IR_typePrimer_OFDM,NIR_typePrimer_OFDM] = PrimerActuator( ...
    M_SNR_dB,...                    % Signal-to-Noise Ratio (dB)
    G_OFDM,...                      % Transmit Matrix
    Q_OFDM,...                      % Receive matrix
    GP_OFDM,...                     % Transmit Matrix at the Pilot
    QP_OFDM,...                     % Receive matrix at the pilot
    R_H,...                         % Channel Correlation Matrix
    N,...                           % Total number of samples
    NrPilotSymbols, ...             % Pilot symbol number
    ZeroThresholdSparse, ...        % Sparse Threshold
    R_hP_est_OFDM, ...              % Correlation matrix for channels estimated at the pilot
    R_hP_est_noInterference_OFDM ...% Correlation matrix for channels estimated at the pilot (no Interference)
    );
% FBMC
[IR_typePrimer,NIR_typePrimer] = PrimerActuator( ...
    M_SNR_dB,...                         % Signal-to-Noise Ratio (dB)
    G_FBMC,...                           % Transmit Matrix
    Q_FBMC,...                           % Receive matrix
    GP_FBMC,...                          % Transmit Matrix at the Pilot
    QP_FBMC,...                          % Receive matrix at the pilot
    R_H,...                              % Channel Correlation Matrix
    N,...                                % Total number of samples
    NrPilotSymbols, ...                  % Pilot symbol number
    ZeroThresholdSparse, ...             % Sparse Threshold
    R_hP_est_FBMC, ...                   % Correlation matrix for channels estimated at the pilot
    R_hP_est_noInterference_FBMC ...     % Correlation matrix for channels estimated at the pilot (no Interference)
    );

%% Start Simulation
disp('Monte Carlo Simulation ...');
Velocity_kmh  = [500 560 600];                         % Velocity in km/h. Note that [mph]*1.6=[kmh] and [m/s]*3.6=[kmh] 
for i_V = 1:length(Velocity_kmh)
    tic
    %% Channel Model Object
    ChannelModel = FastFading(...
    SamplingRate,...                                   % Sampling rate (Samples/s)
    PowerDelayProfile,...                              % Power delay profile, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate)
    N,...                                              % Number of total samples
    Velocity_kmh(i_V)/3.6*CarrierFrequency/2.998e8,... % Maximum Doppler shift: Velocity_kmh/3.6*CarrierFrequency/2.998e8
    'Jakes',...                                        % Which Doppler model: 'Jakes', 'Uniform', 'Discrete-Jakes', 'Discrete-Uniform'. For "Discrete-", we assume a discrete Doppler spectrum to improve the simulation time. This only works accuratly if the number of samples and the velocity is sufficiently large
    200,...                                            % Number of paths for the WSSUS process. Only relevant for a 'Jakes' and 'Uniform' Doppler spectrum
    1,...                                              % Number of transmit antennas
    1,...                                              % Number of receive antennas
    0 ...                                              % Gives a warning if the predefined delay taps of the channel do not fit the sampling rate. This is usually not much of a problem if they are approximatly the same.
    );

for i_rep = 1:NrRepetitions
    %% Update Channel
    ChannelModel.NewRealization;
    %% Binary Data
    BinaryDataStream_OFDM = randi([0 1],NrDataSymbols_OFDM*log2(QAM.ModulationOrder),1);
    BinaryDataStream_FBMC = randi([0 1],GuardSymbolStructure.NrDataSymbols*log2(PAM.ModulationOrder),1);    
    %% Data Symbols
    xD_OFDM = QAM.Bit2Symbol(BinaryDataStream_OFDM);
    xD_FBMC = PAM.Bit2Symbol(BinaryDataStream_FBMC);   
    %% Pilot Symbols
    xP_OFDM = QAM.SymbolMapping(randi(QAM.ModulationOrder,NrPilotSymbols,1));
    xP_OFDM = xP_OFDM./abs(xP_OFDM);

    xP_FBMC = PAM.SymbolMapping(randi(PAM.ModulationOrder,GuardSymbolStructure.NrPilotSymbols,1));
    xP_FBMC = xP_FBMC./abs(xP_FBMC);
    %% Transmitted Data Symbols (Map bin to symbol)
    x_OFDM = PilotMapping_OFDM*[xP_OFDM;xD_OFDM];
    x_FBMC = GuardSymbolStructure.PrecodingMatrix*[xP_FBMC;xD_FBMC];
    x_FBMC = reshape(x_FBMC,FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols);            
    %% Transmitted Signal (time domain)
    s_OFDM  = G_OFDM*x_OFDM(:);
    s_FBMC  = FBMC.Modulation(x_FBMC);  % Same as "s_FBMC = G_FBMC*x_FBMC(:);"  
    %% Channel
    ConvolutionMatrix = ChannelModel.GetConvolutionMatrix{1};
    r_OFDM_noNoise  = ConvolutionMatrix*s_OFDM;
    r_FBMC_noNoise  = ConvolutionMatrix*s_FBMC;
    %% Transmission Matrix
    D_OFDM = Q_OFDM'*ConvolutionMatrix*G_OFDM;
    D_FBMC = Q_FBMC'*ConvolutionMatrix*G_FBMC;    
    %% One-Tap Channel (for perfect channel knowledge)
    h_OFDM = diag(D_OFDM);
    h_FBMC = diag(D_FBMC);
              
    for i_SNR = 1:length(M_SNR_dB)
        %% Add Noise
        SNR_dB  = M_SNR_dB(i_SNR);
        Pn_time = SamplingRate/(F*L)*10^(-SNR_dB/10);
        noise   = sqrt(Pn_time/2)*(randn(size(s_OFDM))+1j*randn(size(s_OFDM)));
        r_OFDM  = r_OFDM_noNoise+noise;
        r_FBMC  = r_FBMC_noNoise+noise;
        %% Demodulate FBMC signal
        y_OFDM  = Q_OFDM'*r_OFDM; % Same as "OFDM.Demodulation(r_OFDM)" 
        y_FBMC  = Q_FBMC'*r_FBMC; % Same as "FBMC.Demodulation(r_FBMC)"
        %% Channel Estimation at Pilot Position
        hP_est_OFDM = y_OFDM(PilotMatrix_OFDM==1)./xP_OFDM/sqrt(Rho_OFDM);
        hP_est_FBMC = y_FBMC(PilotMatrix_FBMC==1)./xP_FBMC/sqrt(Rho_FBMC);
       
        %% Estimate Transmit Matix 
        D_OFDM_est = sum(bsxfun(@times,...
            reshape(full(NIR_typePrimer_OFDM(:,i_SNR)),size(G_OFDM,2),size(G_OFDM,2),NrPilotSymbols),...
            reshape(hP_est_OFDM,1,1,[])),3);     
        h_est_OFDM = diag(D_OFDM_est);
        x_est_OneTapEqualizer_OFDM = y_OFDM./h_est_OFDM;
        xD_est_OneTapEqualizer_OFDM = x_est_OneTapEqualizer_OFDM(PilotMatrix_OFDM(:)==0)./sqrt(DataPowerReduction_OFDM); 

        D_FBMC_est = sum(bsxfun(@times,...
            reshape(full(NIR_typePrimer(:,i_SNR)),size(G_FBMC,2),size(G_FBMC,2),NrPilotSymbols),...
            reshape(hP_est_FBMC,1,1,[])),3);   
        h_est_FBMC = diag(D_FBMC_est);  
        x_est_OneTapEqualizer_FBMC  = y_FBMC./h_est_FBMC;
        xD_est_OneTapEqualizer_FBMC = real(x_est_OneTapEqualizer_FBMC(GuardPilotMatrix_FBMC(:)==0)./sqrt(GuardSymbolStructure.DataPowerReduction));
        
        %% Improved Channel Estimation and Data Detection    
        xD_est_OFDM_Temp     = xD_est_OneTapEqualizer_OFDM;     % initialize with one tap estimates    
        [H_Est_OFDM,h_est_OFDM_Temp,y_OFDM_InterferenceCancellation]  = ImprovedChannelEstimation( ...
                  i_SNR,...                     % Noise Index 
                  NrIterations, ...             % Number of iterations
                  NrPilotSymbols, ...           % Pilot symbol number
                  G_OFDM, ...                   % Transmitted Matrix
                  D_OFDM_est, ...               % The estimated channel matrix
                  y_OFDM, ...                   % The received symbol
                  xD_est_OFDM_Temp,...          % The estimated data symbols
                  xP_OFDM, ...                  % The Transmitted pilot symbols
                  Rho_OFDM, ...               % Pilot-to-data power offset
                  PilotMatrix_OFDM, ...         % Pilot Matrix
                  PilotMapping_OFDM, ...        % Auxiliary matrix
                  NIR_typePrimer_OFDM, ...             % B_type Primer vector
                  IR_typePrimer_OFDM,...              % Interference-free A_type Primer vector
                  QAM ...                       % Modulation Method
                  );

        xD_est_FBMC_Temp = xD_est_OneTapEqualizer_FBMC; % initialize with one tap estimates      
        [H_Est_FBMC,h_est_FBMC_Temp,y_FBMC_InterferenceCancellation]  = ImprovedChannelEstimation( ...
                  i_SNR,...                            % Noise Index
                  NrIterations, ...                    % Number of iterations
                  NrPilotSymbols, ...                  % Pilot symbol number
                  G_FBMC, ...                          % Transmitted Matrix
                  D_FBMC_est, ...                  % The estimated channel matrix
                  y_FBMC(:), ...                   % The received symbol
                  xD_est_FBMC_Temp,...             % The estimated data symbols
                  xP_FBMC, ...                         % The Transmitted pilot symbols
                  Rho_FBMC, ...                       % Pilot-to-data power offset
                  PilotMatrix_FBMC, ...                % Pilot Matrix
                  GuardSymbolStructure.PrecodingMatrix, ... % Auxiliary matrix
                  NIR_typePrimer, ...                 % B_type Primer vector
                  IR_typePrimer,...                  % Interference-free A_type Primer vector
                  PAM ...                              % Modulation Method
                  );
    %% Channel H NMSE 
       NMSE_OFDM_Imp(i_SNR,i_rep) = 10*log10((norm(abs(H_Est_OFDM)-abs(D_OFDM)))^2/(norm(D_OFDM))^2);
       NMSE_FBMC_Imp(i_SNR,i_rep) = 10*log10((norm(abs(H_Est_FBMC)-abs(D_FBMC)))^2/(norm(D_FBMC))^2);
    end
    TimeNeededSoFar = toc;
    disp([int2str(i_rep/NrRepetitions*100) '% Completed! Time Left: ' int2str(TimeNeededSoFar/i_rep*(NrRepetitions-i_rep)/60) 'min, corresponding to approx. '  int2str(TimeNeededSoFar/i_rep*(NrRepetitions-i_rep)/3600) 'hour']);
end
    switch Velocity_kmh(i_V)
        case 500
            NMSE_OFDM_500 = NMSE_OFDM_Imp;
            NMSE_FBMC_500 = NMSE_FBMC_Imp;
        case 560
            NMSE_OFDM_560 = NMSE_OFDM_Imp;
            NMSE_FBMC_560 = NMSE_FBMC_Imp;
        case 600
            NMSE_OFDM_600 = NMSE_OFDM_Imp;
            NMSE_FBMC_600 = NMSE_FBMC_Imp;
    end
    disp(['Complet ' i_V/length(Velocity_kmh)]);
end
%% Plot results
Markersize = 9; 
LineWidth  = 1.3;
figure(1);
subplot(1, 2, 1);
hold off; 
plot(M_SNR_dB, mean(NMSE_OFDM_500,2),      '-o','Color' ,[1 0 0]*0.75,'LineWidth',LineWidth,'Markersize',Markersize);
hold on; grid on;
plot(M_SNR_dB, mean(NMSE_OFDM_560,2),      '-*','Color' ,[0 1 0]*0.75,'LineWidth',LineWidth,'Markersize',Markersize);
plot(M_SNR_dB, mean(NMSE_OFDM_600,2),      '-d','Color' ,[0 0 1]*0.75,'LineWidth',LineWidth,'Markersize',Markersize);
ylim([-25 0]);
title('OFDM, MIR-${\bf R}_{\bf H}$-560.','Interpreter','latex')
legend({'$\hat{\mathcal{H}}$ 500 km/h','$\hat{\mathcal{H}}$ 560 km/h', ...
        '$\hat{\mathcal{H}}$ 600 km/h' },'Interpreter','latex','Location','northeast');
ylabel('NMSE (dB)');
xlabel('Signal-to-Noise Ratio (dB)');
set(gca,'XTick',[0 10 20 30 40],'FontName','Times New Roman','FontSize',15,'GridLineStyle','--','GridLineWidth',1.2);

subplot(1, 2, 2);
hold off; 
plot(M_SNR_dB, mean(NMSE_FBMC_500,2),      '-o','Color' ,[1 0 0]*0.75,'LineWidth',LineWidth,'Markersize',Markersize);
hold on; grid on;
plot(M_SNR_dB, mean(NMSE_FBMC_560,2),      '-*','Color' ,[0 1 0]*0.75,'LineWidth',LineWidth,'Markersize',Markersize);
plot(M_SNR_dB, mean(NMSE_FBMC_600,2),      '-d','Color' ,[0 0 1]*0.75,'LineWidth',LineWidth,'Markersize',Markersize);
ylim([-25 0]);
title('FBMC, MIR-${\bf R}_{\bf H}$-560.','Interpreter','latex')
legend({'$\hat{\mathcal{H}}$ 500 km/h','$\hat{\mathcal{H}}$ 560 km/h', ...
        '$\hat{\mathcal{H}}$ 600 km/h' },'Interpreter','latex','Location','northeast');
ylabel('NMSE (dB)');
xlabel('Signal-to-Noise Ratio (dB)');
set(gca,'XTick',[0 10 20 30 40],'FontName','Times New Roman','FontSize',15,'GridLineStyle','--','GridLineWidth',1.2);

