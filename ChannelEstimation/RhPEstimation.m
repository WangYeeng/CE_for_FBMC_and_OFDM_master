function [R_hP_est,R_hP_est_noInterference] = RhPEstimation( ...
    SNR_dB, ...          % Signal-to-Noise Ratio (dB)
    Bandwidth, ...       % Bandwidth
    Kappa, ...           % Pilot-to-data power offset
    GP_TX, ...           % Transmit Matrix at the Pilot
    GP_RX, ...           % Receive matrix at the pilot
    G_Meth, ...          % Transmit Matrix
    R_H, ...             % Channel Correlation Matrix
    SamplingRate, ...    % Sampling Rate
    NrPilotSymbols, ...  % Pilot symbol number
    NrSamplesTotal ...   % Total number of samples
    )
%Calculate the channel correlation matrix at the pilot position, 
% including the interference-free case.
%% Calculate Correlation Matrices
R_hP  = nan(NrPilotSymbols,NrPilotSymbols);
for j_pilot = 1:NrPilotSymbols
    R_hP(:,j_pilot) = sum((GP_RX'*reshape(R_H*kron(GP_TX(:,j_pilot).',GP_RX(:,j_pilot)')',NrSamplesTotal,NrSamplesTotal)).*(GP_TX.'),2);    
end

R_hP_est_noNoise = R_hP;
for i_pilots = 1: NrPilotSymbols
    Temp = kron(sparse(eye(NrSamplesTotal)),GP_RX(:,i_pilots)')/sqrt(Kappa);
    R_hP_est_noNoise(i_pilots,i_pilots)    =abs(sum(sum((G_Meth.'*(Temp*R_H*Temp')).*G_Meth',2)));           
end 

R_hP_est = repmat(R_hP_est_noNoise,[1 1 length(SNR_dB)]);
for i_SNR = 1:length(SNR_dB)
    SNR = SNR_dB(i_SNR);
    Pn_time = SamplingRate/(Bandwidth)*10^(-SNR/10);   
    for i_pilots = 1: NrPilotSymbols     
        R_hP_est(i_pilots,i_pilots,i_SNR) = R_hP_est_noNoise(i_pilots,i_pilots)+Pn_time*GP_RX(:,i_pilots)'*GP_RX(:,i_pilots)/Kappa;           
    end    
end
R_hP_est_noInterference = R_hP_est-repmat((R_hP_est_noNoise-R_hP),[1 1 length(SNR_dB)]);

end