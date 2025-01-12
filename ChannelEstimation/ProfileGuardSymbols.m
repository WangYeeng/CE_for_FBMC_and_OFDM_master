classdef ProfileGuardSymbols < handle 
   % =====================================================================        
   % This MATLAB class represents an implementation of Profile Guard Symbol
   % for pilot aided channel estimation in FBMC.
   % =====================================================================    
   properties (SetAccess = private)
       Method
       PilotMatrix
       PrecodingMatrix
       PostCodingChannelMatrix
       NrDataSymbols
       NrPilotSymbols
       NrGuardSymbols
       NrTransmittedSymbols
       PilotToDataPowerOffset
       GuardToDataPowerOffset
       DataPowerReduction
       SIR_dB
       ConsideredInterferenceMatrix
   end
   
   
   methods
      % Class constructor, define default values. 
      function obj = ProfileGuardSymbols(varargin)
          % Initialize parameters, set default values
          Method                        = varargin{1};                  % Cancellation method, either 'Coding' or 'Auxiliary'
          PilotMatrix                   = varargin{2};                  % PilotMatrix, 0 = Data, 1 = Pilot, -1 = Auxiliary symbol                
          FBMCMatrix                    = varargin{3};                  % FBMC transmission matrix D, i.e., y = D*x with x transmitted data symbols and y received data symbols (before equalization)
          NrCanceledInterferersPerPilot = varargin{4};                  % Number of neighboring time-frequency positions which are canceled. The higher the number, the lower the SIR but the higher the complexity. For the coding approach, the canceled symbols must not overlap
          PilotToDataPowerOffset        = varargin{5};                  % Pilot to data power offset. 2 guarantees that the SNR is the same at pilot position and at data position => fair comparision.

          
          % Abs Interference Matrix, same as InterferenceMatrix = FBMC.GetInterferenceMatrix
          InterferenceMatrix_11         = abs(reshape(FBMCMatrix(:,1),size(PilotMatrix)));
          InterferenceMatrix_End1       = abs(reshape(FBMCMatrix(:,size(PilotMatrix,1)),size(PilotMatrix)));
          InterferenceMatrix_1End       = abs(reshape(FBMCMatrix(:,numel(PilotMatrix)-size(PilotMatrix,1)+1),size(PilotMatrix)));
          InterferenceMatrix_EndEnd     = abs(reshape(FBMCMatrix(:,numel(PilotMatrix)),size(PilotMatrix)));
          InterferenceMatrix            = [[InterferenceMatrix_EndEnd;InterferenceMatrix_1End(2:end,:)],[InterferenceMatrix_End1(:,2:end);InterferenceMatrix_11(2:end,2:end)]];         
          
          
          switch Method
              case 'Guard'                  
                  NrPilotSymbols        = sum(PilotMatrix(:)==1);
                  NrDataSymbols         =  sum(PilotMatrix(:)==0);
                  NrGuardSymbols        = sum(PilotMatrix(:)==-1);

                  PseudoInvers      = pinv(FBMCMatrix(PilotMatrix(:)==1,PilotMatrix(:)==-1));
                  GuaMatrixPilots   = PseudoInvers*(eye(NrPilotSymbols)-FBMCMatrix(PilotMatrix(:)==1,PilotMatrix(:)==1));
                  GuaMatrixData     = -PseudoInvers*FBMCMatrix(PilotMatrix(:)==1,PilotMatrix(:)==0);


                  GuardMatrix = zeros(numel(PilotMatrix),numel(PilotMatrix)-NrGuardSymbols);
                  GuardMatrix(PilotMatrix==-1,1:NrPilotSymbols)     = GuaMatrixPilots;
                  GuardMatrix(PilotMatrix==-1,NrPilotSymbols+1:end) = GuaMatrixData;
                  GuardMatrix(PilotMatrix==1,1:NrPilotSymbols)      = eye(NrPilotSymbols)*sqrt(PilotToDataPowerOffset);
                  GuardMatrix(PilotMatrix==0,NrPilotSymbols+1:end)  = eye(NrDataSymbols);

                  if NrCanceledInterferersPerPilot>0            
                      [SortedInterferenceValues]    = sort(abs(InterferenceMatrix(:)),'descend');
                      ConsideredInterference_temp   = abs(FBMCMatrix(PilotMatrix(:)==1,:))>=SortedInterferenceValues(NrCanceledInterferersPerPilot+1);
                      Temp_PilotNumber(1,1,:)       = -(1:NrPilotSymbols);
                      ConsideredInterference        = sum(bsxfun(@times,reshape(ConsideredInterference_temp',size(PilotMatrix,1),size(PilotMatrix,2),NrPilotSymbols),Temp_PilotNumber),3);
                      ConsideredInterference(PilotMatrix(:)==1) = 1:NrPilotSymbols; 
                      Index_Pilots  = ConsideredInterference;
                      Index_Data    =  ConsideredInterference;
                      Index_Pilots(PilotMatrix(:)<1)=[];
                      Index_Data(PilotMatrix(:)~=0)=[];

                      GuardMatrix(PilotMatrix(:)==-1,[Index_Pilots,Index_Data]==0)=0;            
                  else
                      ConsideredInterference = 'All';
                  end

                  % Normalize
                  DataPowerReduction    = (numel(PilotMatrix)/sum(sum(GuardMatrix.*conj(GuardMatrix),2),1));
                  GuardMatrix           = GuardMatrix*sqrt(DataPowerReduction);


                  FBMCMatrix_Temp   = FBMCMatrix(PilotMatrix==1,:)*GuardMatrix;
                  SIR_dB= nan(NrPilotSymbols,1);
                  for i_pilot = 1: NrPilotSymbols
                      SIR_dB(i_pilot) = 10*log10(abs(FBMCMatrix_Temp(i_pilot,i_pilot)).^2/(sum(abs(FBMCMatrix_Temp(i_pilot,:)).^2)-abs(FBMCMatrix_Temp(i_pilot,i_pilot)).^2));
                  end

                  Power = diag(GuardMatrix*GuardMatrix');

                  GuardToDataPowerOffset = mean(Power(PilotMatrix(:)==-1))./mean(Power(PilotMatrix(:)==0));
                
                  PrecodingMatrix = GuardMatrix; 
                  obj.PostCodingChannelMatrix = nan;            

              otherwise
                  error('Method must be  ''Guard''!');
          end  
          
          % Set Properties
          obj.Method                        = Method;
          obj.PilotMatrix                   = PilotMatrix;
          obj.PrecodingMatrix               = PrecodingMatrix;
          obj.NrDataSymbols                 = NrDataSymbols;
          obj.NrPilotSymbols                = NrPilotSymbols; 
          obj.NrGuardSymbols                = NrGuardSymbols;
          obj.NrTransmittedSymbols          = size(PrecodingMatrix,1);
          obj.PilotToDataPowerOffset        = PilotToDataPowerOffset;
          obj.GuardToDataPowerOffset        = GuardToDataPowerOffset;
          obj.DataPowerReduction            = DataPowerReduction;
          obj.SIR_dB                        = SIR_dB;
          obj.ConsideredInterferenceMatrix  = ConsideredInterference;     
          
      end
           
   end 
end

