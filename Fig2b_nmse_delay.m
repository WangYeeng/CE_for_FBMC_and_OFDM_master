% =========================================================================   
% (c) 2024 Ying Wang, wangyingstu@163.com
% =========================================================================  
clc; clear; close all;
addpath('./MultiCarrierSystem', ...
        './ChannelEstimation','./Channel');
load("R_H_560.mat");
%% Parameters
% Simulation
M_SNR_dB                  = 25;                      % Signal-to-Noise Ratio in dB
NrRepetitions             = 1000;                    % Number of Monte Carlo repetition (different channel realizations)             
ZeroThresholdSparse       = 10;                      % Set some matrix values, smaller than "10^(-ZeroThresholdSparse)", to zero.

% FBMC and OFDM parameters
L                         = 24;                     % Number of Subcarriers, one resource block consists of 12 subcarriers (and 0.5ms in time)
F                         = 15e3;                   % Subcarrier Spacing in Hz
SamplingRate              = F*L*5;                  % Sampling rate in (Samples/s)
CarrierFrequency          = 2.5e9;                  % Carrier Frequency
Cycliclength              = 1/15e3/14;
NrSubframes               = 1;                      % Number of subframes. One subframe requires 1ms for F=15kHz.                             
QAM_ModulationOrder       = 16;                     % QAM signal constellation order, 4, 16, 64, 256, 1024,...

% Channel estimation parameters
NrIterations              = 5;                      % Number of iterations for interference cancellation scheme.
PilotToDataPowerOffset    = 2.15;                  
% Channel
Velocity_kmh              = 500;                    % Velocity in km/h. Note that [mph]*1.6=[kmh] and [m/s]*3.6=[kmh]        
Delay                     = [50 150 250 350 450 550];
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

N = OFDM.Nr.SamplesTotal;
%% QAM Object
QAM = SignalConstellation(QAM_ModulationOrder,'QAM');

%% Pilot Matrices, 0="Data", 1=Pilot
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
%%=---------------------------Same as-------------------------------------- 
% Rowinterval = 6;
% Colinterval = 3.5;
% 
% PilotMatrix_OFDM = zeros(OFDM.Nr.Subcarriers,14);
% PilotMatrix_OFDM(2:2*Rowinterval:end,2:2*Colinterval:end)=1;
% PilotMatrix_OFDM(5:2*Rowinterval:end,6:2*Colinterval:end)=1;
% PilotMatrix_OFDM(8:2*Rowinterval:end,2+1:2*Colinterval:end)=1;
% PilotMatrix_OFDM(11:2*Rowinterval:end,6+1:2*Colinterval:end)=1;
%%-------------------------------------------------------------------------
PilotMatrix_OFDM = repmat(PilotMatrix_OFDM,[1 NrSubframes]);

%% Cancel Imaginary Interference At Pilot Position Object (Precoding Matrix)
NrPilotSymbols     = sum(PilotMatrix_OFDM(:)==1);
NrDataSymbols_OFDM = sum(PilotMatrix_OFDM(:)==0);

PilotMapping_OFDM = zeros(numel(PilotMatrix_OFDM),numel(PilotMatrix_OFDM));
PilotMapping_OFDM(PilotMatrix_OFDM(:)==1,1:NrPilotSymbols) = sqrt(PilotToDataPowerOffset)*eye(NrPilotSymbols);
PilotMapping_OFDM(PilotMatrix_OFDM(:)==0,NrPilotSymbols+1:end) = eye(NrDataSymbols_OFDM);
PilotMapping_OFDM = PilotMapping_OFDM/sqrt(mean(diag(PilotMapping_OFDM*PilotMapping_OFDM')));
DataPowerReduction_OFDM = numel(PilotMatrix_OFDM)/(NrPilotSymbols*PilotToDataPowerOffset+NrDataSymbols_OFDM);

Rho_OFDM     = PilotToDataPowerOffset*DataPowerReduction_OFDM;
Overhead_OFDM  = (NrPilotSymbols)/(OFDM.Nr.Subcarriers*OFDM.Nr.MCSymbols) + (Cycliclength)/OFDM.PHY.TimeSpacing;

%% Precalculate Transmit and Receive Matrices
G_OFDM = OFDM.GetTXMatrix;
Q_OFDM = (OFDM.GetRXMatrix)';
GP_OFDM = G_OFDM(:,PilotMatrix_OFDM(:)==1);
QP_OFDM = Q_OFDM(:,PilotMatrix_OFDM(:)==1);
G_OFDM_PilotMapping = G_OFDM*PilotMapping_OFDM;

disp('Calculate Correlation Matrixe...');
%% Calculate Correlation Matrix of Pilot Estimates for OFDM
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

disp('Calculate induction factor Solution and that of no Interference...');
%% Calculate MMSE Estimation Matrix and that of no Interference for OFDM
[IR_typePrimer,NIR_typePrimer] = PrimerActuator( ...
    M_SNR_dB,...                    % Signal-to-Noise Ratio (dB)
    G_OFDM,...                      % Transmit Matrix
    Q_OFDM,...                      % Receive matrix
    GP_OFDM,...                     % Transmit Matrix at the Pilot
    QP_OFDM,...                     % Receive matrix at the pilot
    R_H,...                      % Channel Correlation Matrix
    N,...                           % Total number of samples
    NrPilotSymbols, ...             % Pilot symbol number
    ZeroThresholdSparse, ...        % Sparse Threshold
    R_hP_est_OFDM, ...              % Correlation matrix for channels estimated at the pilot
    R_hP_est_noInterference_OFDM ...% Correlation matrix for channels estimated at the pilot (no Interference)
    );

%% Start Simulation
tic
disp('Monte Carlo Simulation ...');
for i_rep = 1:NrRepetitions
    %% Binary Data
    BinaryDataStream_OFDM     = randi([0 1],NrDataSymbols_OFDM*log2(QAM.ModulationOrder),1);
       
    %% Data Symbols
    xD_OFDM     = QAM.Bit2Symbol(BinaryDataStream_OFDM);
       
    %% Pilot Symbols
    xP_OFDM = QAM.SymbolMapping(randi(QAM.ModulationOrder,NrPilotSymbols,1));
    xP_OFDM = xP_OFDM./abs(xP_OFDM);
    
    %% Transmitted Data Symbols (Map bin to symbol)
    x_OFDM     = PilotMapping_OFDM*[xP_OFDM;xD_OFDM];
                
    %% Transmitted Signal (time domain)
    s_OFDM     = G_OFDM*x_OFDM(:);
for i_delay = 1:length(Delay)       
    %% Channel
    %% Channel Model Object
    variableName = ['TDL-A_', num2str(Delay(i_delay)), 'ns'];
    ChannelModel = FastFading(...
             SamplingRate,...                                   % Sampling rate (Samples/s)
             num2str(variableName),...                              % Power delay profile, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate)
             N,...                                              % Number of total samples
             Velocity_kmh/3.6*CarrierFrequency/2.998e8,...      % Maximum Doppler shift: Velocity_kmh/3.6*CarrierFrequency/2.998e8
             'Jakes',...                                        % Which Doppler model: 'Jakes', 'Uniform', 'Discrete-Jakes', 'Discrete-Uniform'. For "Discrete-", we assume a discrete Doppler spectrum to improve the simulation time. This only works accuratly if the number of samples and the velocity is sufficiently large
             200,...                                            % Number of paths for the WSSUS process. Only relevant for a 'Jakes' and 'Uniform' Doppler spectrum
             1,...                                              % Number of transmit antennas
             1,...                                              % Number of receive antennas
             0 ...                                              % Gives a warning if the predefined delay taps of the channel do not fit the sampling rate. This is usually not much of a problem if they are approximatly the same.
             );

    ConvolutionMatrix = ChannelModel.GetConvolutionMatrix{1};
    r_OFDM_noNoise     = ConvolutionMatrix*s_OFDM;
    
    %% Transmission Matrix
    D_OFDM = Q_OFDM'*ConvolutionMatrix*G_OFDM;
        
    %% One-Tap Channel (for perfect channel knowledge)
    h_OFDM = diag(D_OFDM);
    h_OFDM_perfect  = reshape( sum((Q_OFDM'*ConvolutionMatrix).*G_OFDM.',2), L, [] );
           
    for i_SNR = 1:length(M_SNR_dB)
        %% Add Noise
        SNR_dB  = M_SNR_dB(i_SNR);
        Pn_time = SamplingRate/(F*L)*10^(-SNR_dB/10);
        noise   = sqrt(Pn_time/2)*(randn(size(s_OFDM))+1j*randn(size(s_OFDM)));
        r_OFDM     = r_OFDM_noNoise+noise;

        %% Demodulate FBMC signal
        y_OFDM             = Q_OFDM'*r_OFDM; % Same as "OFDM.Demodulation(r_OFDM)" 

        %% Channel Estimation at Pilot Position
        hP_est_OFDM     = y_OFDM(PilotMatrix_OFDM==1)./xP_OFDM/sqrt(Rho_OFDM);
 
        %% Estimate Transmit Matix 
        D_OFDM_est = sum(bsxfun(@times,...
            reshape(full(NIR_typePrimer(:,i_SNR)),size(G_OFDM,2),size(G_OFDM,2),NrPilotSymbols),...
            reshape(hP_est_OFDM,1,1,[])),3);     
        h_est_OFDM = diag(D_OFDM_est);
        x_est_OneTapEqualizer_OFDM = y_OFDM./h_est_OFDM;
        xD_est_OneTapEqualizer_OFDM = x_est_OneTapEqualizer_OFDM(PilotMatrix_OFDM(:)==0)./sqrt(DataPowerReduction_OFDM);   
        h_OFDM_Inter    = ChannelEstimation_OFDM.ChannelInterpolation(hP_est_OFDM);
 
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
                  NIR_typePrimer, ...             % B_type Primer vector
                  IR_typePrimer,...              % Interference-free A_type Primer vector
                  QAM ...                       % Modulation Method
                  );
    %% Channel H NMSE 
       NMSE_OFDM_LS(i_delay,i_rep)  = 10*log10((norm(abs(D_OFDM_est)-abs(D_OFDM)))^2/(norm(D_OFDM))^2);
       NMSE_OFDM_Imp(i_delay,i_rep) = 10*log10((norm(abs(H_Est_OFDM)-abs(D_OFDM)))^2/(norm(D_OFDM))^2);
       %% Channel_Temp h NMSE
       NMSE_OFDM_LS_Temp(i_delay,i_rep)    = 10*log10((norm(abs(h_est_OFDM)-abs(h_OFDM)))^2/(norm(h_OFDM))^2);
       NMSE_OFDM_Imp_Temp(i_delay,i_rep)   = 10*log10((norm(abs(h_est_OFDM_Temp)-abs(h_OFDM)))^2/(norm(h_OFDM))^2);
       NMSE_OFDM_Inter_Temp(i_delay,i_rep) = 10*log10((norm(abs(h_OFDM_Inter(:))-abs(h_OFDM)))^2/(norm(h_OFDM))^2);
       
    end
end
    TimeNeededSoFar = toc;
    disp([int2str(i_rep/NrRepetitions*100) '% Completed! Time Left: ' int2str(TimeNeededSoFar/i_rep*(NrRepetitions-i_rep)/60) 'min, corresponding to approx. '  int2str(TimeNeededSoFar/i_rep*(NrRepetitions-i_rep)/3600) 'hour']);
end

%% Plot results
Markersize = 10; 
LineWidth  = 1.3;
figure(1);
hold off; 
plot(Delay, mean(NMSE_OFDM_Inter_Temp,2),'--','Color' ,[1 0 0]*0.75,'LineWidth',LineWidth,'Markersize',Markersize);
hold on; grid on;box on;
plot(Delay, mean(NMSE_OFDM_LS,2),       '- *','Color' ,[1 0 1]*0.75,'LineWidth',LineWidth,'Markersize',Markersize);
plot(Delay, mean(NMSE_OFDM_Imp,2),      '- o','Color' ,[0 1 0]*0.70,'LineWidth',LineWidth,'Markersize',Markersize);
plot(Delay, mean(NMSE_OFDM_LS_Temp,2),  '- d','Color' ,[0 0 1]*0.80,'LineWidth',LineWidth,'Markersize',Markersize);
plot(Delay, mean(NMSE_OFDM_Imp_Temp,2), '-. ' ,'Color' ,[1 1 1]*0.10,'LineWidth',LineWidth,'Markersize',Markersize);
ylabel('Normalized Mean Square Error (dB)');
xlabel('RMS Delay (ns)');
ylim([-35 -10]);
legend({'$\hat{h}$ est. by LS','$\hat{\mathcal{H}}$est. by NIR','$\hat{\mathcal{H}}$est. by IR', ...
        '$\hat{h}$ est. by NIR','$\hat{h}$ est. by IR'},'Interpreter','latex','Location','southeast','NumColumns', 3);
set(gca,'Ytick',[-35 -30 -25 -20 -15 -10],'Xtick',[50 150 250 350 450 550],'FontName','Times New Roman','FontSize',15,'GridLineStyle','-.','GridLineWidth',1.2,'LooseInset', [0,0,0,0]);

