function [A_typePrimer,B_typePrimer] = PrimerActuator( ...
    SNR_dB,...                  % Signal-to-Noise Ratio (dB)
    G_TX,...                    % Transmit Matrix
    G_RX,...                    % Receive matrix
    GP_TX,...                   % Transmit Matrix at the Pilot
    GP_RX,...                   % Receive matrix at the pilot
    R_H,...                     % Channel Correlation Matrix
    NrSamplesTotal,...          % Total number of samples
    NrPilotSymbols, ...         % Pilot symbol number
    ZeroThresholdSparse, ...    % Sparse Threshold
    R_hP_est, ...               % Correlation matrix for channels estimated at the pilot         
    R_hP_est_noInterference ... % Correlation matrix for channels estimated at the pilot (no Interference)
    )
% Calculate the PullFactor vector, 
% including the case of no interference.
%% CalculatePullFactor Vector
R_Dij_hP = sparse(size(G_TX,2)^2,NrPilotSymbols);
for i_pilot = 1: NrPilotSymbols
    R_Dij_hP_Temp = reshape(G_RX'*reshape(R_H*kron(GP_TX(:,i_pilot).',GP_RX(:,i_pilot)')',NrSamplesTotal,NrSamplesTotal)*G_TX,[],1);
    R_Dij_hP_Temp(abs(R_Dij_hP_Temp)<10^(-ZeroThresholdSparse))=0;
    R_Dij_hP(:,i_pilot) = R_Dij_hP_Temp;
end

A_typePrimer = sparse(size(G_TX,2)*size(G_TX,2)*NrPilotSymbols,length(SNR_dB));
for i_SNR = 1:length(SNR_dB)
    A_typePrimer_Temp = reshape(R_Dij_hP*pinv(R_hP_est_noInterference(:,:,i_SNR)),size(G_TX,2)*size(G_TX,2)*NrPilotSymbols,1);
    A_typePrimer_Temp(abs(A_typePrimer_Temp)<10^(-ZeroThresholdSparse))=0;  
    A_typePrimer(:,i_SNR) = A_typePrimer_Temp;   
end

B_typePrimer = sparse(size(G_TX,2)*size(G_TX,2)*NrPilotSymbols,length(SNR_dB));
for i_SNR = 1:length(SNR_dB)
    B_typePrimer_Temp = reshape(R_Dij_hP*pinv(R_hP_est(:,:,i_SNR)),size(G_TX,2)*size(G_TX,2)*NrPilotSymbols,1);
    B_typePrimer_Temp(abs(B_typePrimer_Temp)<10^(-ZeroThresholdSparse))=0;
    B_typePrimer(:,i_SNR) = B_typePrimer_Temp;
end
clear R_H;
end

