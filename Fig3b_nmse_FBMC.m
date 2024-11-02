% =========================================================================   
% (c) 2024 Ying Wang, wangyingstu@163.com
% =========================================================================  
clc; clear; close all;
addpath('./MultiCarrierSystem', ...
        './ChannelEstimation', './Channel');
load("R_H_560.mat");
NMSE_For_IIU  = [0.31432823, 0.25118864, 0.16943282, 0.10511886, 0.05543284, 0.02711876, 0.0119433, 0.0053121];
% The above data is from article D. Chen, R. Wang, and T. Jiang, 
% "Channel estimation and pilot symbol optimization based on 
% intrinsic interference utilization for OQAM/FBMC systems," 
% IEEE Transactions on Signal Processing, vol. 69, pp. 4595-4606, 2021.
%% Parameters
% Simulation
M_SNR_dB                  = 0:5:35;                 % Signal-to-Noise Ratio in dB
NrRepetitions             = 1000;                     % Number of Monte Carlo repetition (different channel realizations)             
ZeroThresholdSparse       = 10;                     % Set some matrix values, smaller than "10^(-ZeroThresholdSparse)", to zero.

% FBMC and OFDM parameters
L                         = 24;                     % Number of Subcarriers, one resource block consists of 12 subcarriers (and 0.5ms in time)
K                         = 30;
F                         = 15e3;                   % Subcarrier Spacing in Hz
SamplingRate              = F*L*5;                  % Sampling rate in (Samples/s)
CarrierFrequency          = 2.5e9;                  % Carrier Frequency
NrSubframes               = 1;                      % Number of subframes. One subframe requires 1ms for F=15kHz.                             
QAM_ModulationOrder       = 16;                     % QAM signal constellation order, 4, 16, 64, 256, 1024,...

% Channel estimation parameters
NrIterations              = 5;                      % Number of iterations for interference cancellation scheme.                  
PilotToDataPowerOffset    = 4.3; 
% Channel
Velocity_kmh              = 500;                    % Velocity in km/h. Note that [mph]*1.6=[kmh] and [m/s]*3.6=[kmh]        
PowerDelayProfile         = 'ExtendedVehicularA';   % Channel model: 'VehicularA', 'VehicularB', 'ExtendedVehicularA'

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

%% Check Number of Samples
N = FBMC.Nr.SamplesTotal;
%% PAM and QAM Object
PAM = SignalConstellation(sqrt(QAM_ModulationOrder),'PAM');
QAM = SignalConstellation(   QAM_ModulationOrder   ,'QAM');
%% Pilot Matrices, 0="Data", 1=Pilot
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
%%------------------------------Same as------------------------------------
% Rowinterval = 12;
% Colinterval = 16;
% PilotMatrix_FBMC =  zeros(FBMC.Nr.Subcarriers,K);
% PilotMatrix_FBMC(2:Rowinterval:end,3:Colinterval:end) = 1;
% PilotMatrix_FBMC(5:Rowinterval:end,11:Colinterval:end) = 1;
% PilotMatrix_FBMC(8:Rowinterval:end,3+1:Colinterval:end) = 1;
% PilotMatrix_FBMC(11:Rowinterval:end,11+1:Colinterval:end) = 1;
% PilotMatrix_FBMC = repmat(PilotMatrix_FBMC,[1 NrSubframes]);
%%-------------------------------------------------------------------------
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
    PilotToDataPowerOffset ...                          % Pilot to data power offset
    );

NrPilotSymbols     = sum(PilotMatrix_FBMC(:)==1);
NrDataSymbols_FBMC = sum(PilotMatrix_FBMC(:)==0);
NrGuardSymbols     = GuardSymbolStructure.NrGuardSymbols;
Overhead_FBMC      = (NrPilotSymbols+NrGuardSymbols)/(FBMC.Nr.Subcarriers*FBMC.Nr.MCSymbols);

Rho_FBMC  = GuardSymbolStructure.PilotToDataPowerOffset*GuardSymbolStructure.DataPowerReduction;
%% Channel Model Object
ChannelModel = FastFading(...
    SamplingRate,...                                   % Sampling rate (Samples/s)
    PowerDelayProfile,...                              % Power delay profile, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate)
    N,...                                              % Number of total samples
    Velocity_kmh/3.6*CarrierFrequency/2.998e8,...      % Maximum Doppler shift: Velocity_kmh/3.6*CarrierFrequency/2.998e8
    'Jakes',...                                        % Which Doppler model: 'Jakes', 'Uniform', 'Discrete-Jakes', 'Discrete-Uniform'. For "Discrete-", we assume a discrete Doppler spectrum to improve the simulation time. This only works accuratly if the number of samples and the velocity is sufficiently large
    200,...                                            % Number of paths for the WSSUS process. Only relevant for a 'Jakes' and 'Uniform' Doppler spectrum
    1,...                                              % Number of transmit antennas
    1,...                                              % Number of receive antennas
    1 ...                                              % Gives a warning if the predefined delay taps of the channel do not fit the sampling rate. This is usually not much of a problem if they are approximatly the same.
    );
%% Precalculate Transmit and Receive Matrices
G_FBMC  = FBMC.GetTXMatrix;
Q_FBMC  = (FBMC.GetRXMatrix)';
GP_FBMC = G_FBMC(:,PilotMatrix_FBMC(:)==1);
QP_FBMC = Q_FBMC(:,PilotMatrix_FBMC(:)==1);
G_Gua   = G_FBMC*GuardSymbolStructure.PrecodingMatrix;

disp('Calculate Correlation Matrixe...');
%% Calculate Correlation Matrix of Pilot Estimates for FBMC_Aux
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
%% Calculate MMSE Estimation Matrix and that of no Interference for FBMC_Aux
[IR_typePrimer,NIR_typePrimer] = PrimerActuator( ...
    M_SNR_dB,...                         % Signal-to-Noise Ratio (dB)
    G_FBMC,...                           % Transmit Matrix
    Q_FBMC,...                           % Receive matrix
    GP_FBMC,...                          % Transmit Matrix at the Pilot
    QP_FBMC,...                          % Receive matrix at the pilot
    R_H,...                           % Channel Correlation Matrix
    N,...                                % Total number of samples
    NrPilotSymbols, ...                  % Pilot symbol number
    ZeroThresholdSparse, ...             % Sparse Threshold
    R_hP_est_FBMC, ...               % Correlation matrix for channels estimated at the pilot
    R_hP_est_noInterference_FBMC ... % Correlation matrix for channels estimated at the pilot (no Interference)
    );

%% Start Simulation
tic
disp('Monte Carlo Simulation ...');
for i_rep = 1:NrRepetitions
    %% Update Channel
    ChannelModel.NewRealization;
    
    %% Binary Data
    BinaryDataStream_FBMC = randi([0 1],GuardSymbolStructure.NrDataSymbols*log2(PAM.ModulationOrder),1); 
    %% Data Symbols
    xD_FBMC = PAM.Bit2Symbol(BinaryDataStream_FBMC);
    %% Pilot Symbols
    xP_FBMC = PAM.SymbolMapping(randi(PAM.ModulationOrder,GuardSymbolStructure.NrPilotSymbols,1));
    xP_FBMC = xP_FBMC./abs(xP_FBMC);
    
    %% Transmitted Data Symbols (Map bin to symbol)
    x_FBMC = GuardSymbolStructure.PrecodingMatrix*[xP_FBMC;xD_FBMC];
    x_FBMC = reshape(x_FBMC,FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols);     
    %% Transmitted Signal (time domain)
    % s_FBMC = G_FBMC*x_FBMC(:); 
    s_FBMC = FBMC.Modulation(x_FBMC);% Same as "s_FBMC = G_FBMC*x_FBMC(:);" 
    %% Channel
    ConvolutionMatrix  = ChannelModel.GetConvolutionMatrix{1}; 
    r_FBMC_noNoise = ConvolutionMatrix*s_FBMC;
    %% Transmission Matrix
    D_FBMC = Q_FBMC'*ConvolutionMatrix*G_FBMC;
    %% One-Tap Channel (for perfect channel knowledge)
    h_FBMC = diag(D_FBMC);
    for i_SNR = 1:length(M_SNR_dB)
        %% Add Noise
        SNR_dB  = M_SNR_dB(i_SNR);
        Pn_time = SamplingRate/(F*L)*10^(-SNR_dB/10);
        noise   = sqrt(Pn_time/2)*(randn(size(s_FBMC))+1j*randn(size(s_FBMC)));

        r_FBMC = r_FBMC_noNoise+noise;       
        %% Demodulate FBMC signal
        y_FBMC         = Q_FBMC'*r_FBMC; % Same as "FBMC.Demodulation(r_FBMC_Aux)" 
        %% Channel Estimation at Pilot Position
        hP_est_FBMC = y_FBMC(PilotMatrix_FBMC==1)./xP_FBMC/sqrt(Rho_FBMC);
        h_FBMC_Inter    = ChannelEstimation_FBMC.ChannelInterpolation(hP_est_FBMC);
        %% Estimate Transmit Matix
        D_FBMC_est = sum(bsxfun(@times,...
            reshape(full(NIR_typePrimer(:,i_SNR)),size(G_FBMC,2),size(G_FBMC,2),NrPilotSymbols),...
            reshape(hP_est_FBMC,1,1,[])),3);   
        h_est_FBMC = diag(D_FBMC_est);  

        x_est_OneTapEqualizer_FBMC  = y_FBMC./h_est_FBMC;
        xD_est_OneTapEqualizer_FBMC = real(x_est_OneTapEqualizer_FBMC(GuardPilotMatrix_FBMC(:)==0)./sqrt(GuardSymbolStructure.DataPowerReduction));
        
        %% Improved Channel Estimation and Data Detection
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
       NMSE_FBMC_LS(i_SNR,i_rep)  = 10*log10((norm(abs(D_FBMC_est)-abs(D_FBMC)))^2/(norm(D_FBMC))^2);
       NMSE_FBMC_Imp(i_SNR,i_rep) = 10*log10((norm(abs(H_Est_FBMC)-abs(D_FBMC)))^2/(norm(D_FBMC))^2);
       %% Channel_Temp h NMSE
       NMSE_FBMC_LS_Temp(i_SNR,i_rep)  = 10*log10((norm(abs(h_est_FBMC)-abs(h_FBMC)))^2/(norm(h_FBMC))^2);
       NMSE_FBMC_Imp_Temp(i_SNR,i_rep) = 10*log10((norm(abs(h_est_FBMC_Temp)-abs(h_FBMC)))^2/(norm(h_FBMC))^2);
       NMSE_FBMC_Inter_Temp(i_SNR,i_rep) = 10*log10((norm(abs(h_FBMC_Inter(:))-abs(h_FBMC)))^2/(norm(h_FBMC))^2);
    end
    TimeNeededSoFar = toc;
    disp([int2str(i_rep/NrRepetitions*100) '% Completed! Time Left: ' int2str(TimeNeededSoFar/i_rep*(NrRepetitions-i_rep)/60) 'min, corresponding to approx. '  int2str(TimeNeededSoFar/i_rep*(NrRepetitions-i_rep)/3600) 'hour']);

end

%% Plot results 
Markersize = 10; 
LineWidth  = 1.3;
figure(1); 
hold off; 
plot(M_SNR_dB, 10*log10(NMSE_For_IIU),      '-+','Color' ,[1 1 0]*0.65,'LineWidth',LineWidth,'Markersize',Markersize);
hold on; grid on;
plot(M_SNR_dB, mean(NMSE_FBMC_Inter_Temp,2),'--','Color' ,[1 0 0]*0.75,'LineWidth',LineWidth,'Markersize',Markersize);
plot(M_SNR_dB, mean(NMSE_FBMC_LS,2),       '- *','Color' ,[1 0 1]*0.75,'LineWidth',LineWidth,'Markersize',Markersize);
plot(M_SNR_dB, mean(NMSE_FBMC_Imp,2),      '- o','Color' ,[0 1 0]*0.70,'LineWidth',LineWidth,'Markersize',Markersize);
plot(M_SNR_dB, mean(NMSE_FBMC_LS_Temp,2),  '- d','Color' ,[0 0 1]*0.80,'LineWidth',LineWidth,'Markersize',Markersize);
plot(M_SNR_dB, mean(NMSE_FBMC_Imp_Temp,2), '- ' ,'Color' ,[1 1 1]*0.10,'LineWidth',LineWidth,'Markersize',Markersize);
title('FBMC, 64-QAM, EVA, 500km/h@2.5GHz')
legend({'IIU in [7]','LS Linear interp $\hat{h}$','NIR extrap est. $\hat{\mathcal{H}}$','IR extrap est. $\hat{\mathcal{H}}$', ...
        'NIR extrap est. $\hat{h}$','IR extrap est. $\hat{h}$'},'Interpreter','latex','Location','southwest');
ylabel('Normalized Mean Square Error (dB)');
xlabel('Signal-to-Noise Ratio (dB)');
set(gca,'FontName','Times New Roman','FontSize',12,'GridLineStyle','--','GridLineWidth',1.2,'LooseInset', [0,0,0,0]);
