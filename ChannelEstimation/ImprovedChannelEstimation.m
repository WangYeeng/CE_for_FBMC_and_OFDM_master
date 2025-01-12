function [H_Est_Temp,h_est_imp,y_InterferenceCancellation] = ImprovedChannelEstimation( ...
    i_SNR,...                % Noise Index...
    NrIterations, ...        % Number of iterations
    NrPilotSymbols, ...      % Pilot symbol number 
    G_TX, ...                % Transmitted Matrix
    H_Est, ...               % The estimated channel matrix
    y, ...                   % The received symbol
    xD_est,...               % The estimated data symbols
    xP, ...                  % The Transmitted pilot symbols
    Kappa, ...               % Pilot-to-data power offset
    PilotMatrix, ...         % Pilot Matrix
    AuxiliaryMatrix, ...     % Auxiliary matrix
    B_typePrimer, ...        % B_type Primer vector
    A_typePrimer,...         % Interference-free A_type Primer vector
    ModulationMethod ...     % Modulation Method
    )
        H_Est_Temp  = H_Est;
        h_est_imp   = diag(H_Est);
        for i_iteration = 1:NrIterations
            y_InterferenceCancellation = (y - (H_Est_Temp-diag(h_est_imp))*AuxiliaryMatrix*[xP;ModulationMethod.SymbolQuantization(xD_est)]); 
            % New Channel Estimates at Pilot Positions            
            hP_est_Temp = y_InterferenceCancellation(PilotMatrix==1)./xP/sqrt(Kappa);
            %% Improved Channel Estimation
            H_Est_Temp = EstimateTransmitMatix( ...
                G_TX, ...                         % Transmitted Matrix
                i_SNR,...                         % Noise Index
                i_iteration, ...                  % Interference cancellation iterative Indexing
                NrIterations, ...                 % Number of iterations for interference cancellation
                NrPilotSymbols, ...               % Pilot symbol number 
                hP_est_Temp, ...                  % Channel coefficients at pilots
                B_typePrimer, ...                 % B_type Primer vector
                A_typePrimer ...                  % Interference-free A_type Primer vector
                );
            % One-Tap Channel
            h_est_imp     = diag(H_Est_Temp);
        end
end 

function H_est_Temp = EstimateTransmitMatix( ...
    G_TX, ...            % Transmitted Matrix
    i_SNR,...            % Noise Index
    i_iteration, ...     % Interference cancellation iterative Indexing
    NrIterations, ...    % Number of iterations for interference cancellation
    NrPilotSymbols, ...  % Pilot symbol number
    hP_est_Temp, ...     % Channel coefficients at pilots
    B_typePrimer, ...    % B_type Primer vector
    A_typePrimer ...     % Interference-free A_type Primer vector
    )
% Completion of all channel coefficients with the 
% estimated channel coefficients at the pilots.

  if i_iteration<=NrIterations/2
     H_est_Temp = sum(bsxfun(@times, ...
         reshape(full(B_typePrimer(:,i_SNR)),size(G_TX,2),size(G_TX,2),NrPilotSymbols), ...
         reshape(hP_est_Temp,1,1,[])),3);        
  else
     H_est_Temp = sum(bsxfun(@times, ...
         reshape(full(A_typePrimer(:,i_SNR)),size(G_TX,2),size(G_TX,2),NrPilotSymbols), ...
         reshape(hP_est_Temp,1,1,[])),3);
  end
end