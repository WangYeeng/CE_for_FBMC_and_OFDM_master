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