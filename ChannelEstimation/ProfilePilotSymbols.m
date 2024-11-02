classdef ProfilePilotSymbols < handle 
   % =====================================================================        
   % This MATLAB class represents an implementation of pilot-aided channel
   % estimation, that is, it allows to estimate the channel by
   % interpolation/extrapolation.
   % The following interpolation methods are supported:
   %    1) 'linear','nearest','natural' based on the MATLAB built-in
   %        function "scatteredInterpolant"
   %    2) 'FullAverage' averages over all pilots => assumes a double flat
   %        channel.
   % The pilot pattern can be 'Diamond', 'Rectangular'
   % =====================================================================    

   properties (SetAccess = private)
       NrPilotSymbols
       PilotPattern
       PilotSpacingFrequency
       PilotSpacingTime
       InterpolationMethod
       Implementation
       InterpolationProperties
       PilotMatrix
   end
   
   methods
      % Class constructor, define default values. 
      function obj = ProfilePilotSymbols(varargin)
          % Initialize parameters, set default values  
          obj.PilotPattern = varargin{1};                                   
          obj.InterpolationMethod = varargin{3};
 
          % Generate pilot matrix according to the specified pilot pattern.
          % A zero corresponse to a data symbol, a one to a pilot symbol
          switch obj.PilotPattern
              case 'Rectangular'
                  NrSubcarriers = varargin{2}(1,1);
                  obj.PilotSpacingFrequency = varargin{2}(1,2);
                  NrMCSymbols = varargin{2}(2,1);
                  obj.PilotSpacingTime = varargin{2}(2,2);
                  
                  obj.PilotMatrix = zeros(NrSubcarriers,NrMCSymbols);     
                  obj.PilotMatrix(round(mod(NrSubcarriers-1,obj.PilotSpacingFrequency)/2)+1:obj.PilotSpacingFrequency:NrSubcarriers,round(round(mod(NrMCSymbols-1,obj.PilotSpacingTime)/2)+1:obj.PilotSpacingTime:NrMCSymbols)) = true;             
              case 'Diamond'
                  NrSubcarriers = varargin{2}(1,1);
                  obj.PilotSpacingFrequency = varargin{2}(1,2);
                  NrMCSymbols = varargin{2}(2,1);
                  obj.PilotSpacingTime = varargin{2}(2,2);
                  
                  obj.PilotMatrix = zeros(NrSubcarriers,NrMCSymbols);    
                  
                  FrequencyPositionShift = floor((NrSubcarriers-max([(1:2*obj.PilotSpacingFrequency:NrSubcarriers),(1+1/2*obj.PilotSpacingFrequency:2*obj.PilotSpacingFrequency:NrSubcarriers),(1+obj.PilotSpacingFrequency:2*obj.PilotSpacingFrequency:NrSubcarriers),(1+3/2*obj.PilotSpacingFrequency:2*obj.PilotSpacingFrequency:NrSubcarriers)]))/2)+1;
                  TimePositionShift = floor((NrMCSymbols-max([(1:2*obj.PilotSpacingTime:NrMCSymbols),(1+obj.PilotSpacingTime):2*obj.PilotSpacingTime:NrMCSymbols]))/2)+1;
                  obj.PilotMatrix(FrequencyPositionShift                                     :2*obj.PilotSpacingFrequency:NrSubcarriers,TimePositionShift                           :2*obj.PilotSpacingTime:NrMCSymbols) = 1;
                  obj.PilotMatrix(FrequencyPositionShift+round(1/2*obj.PilotSpacingFrequency):2*obj.PilotSpacingFrequency:NrSubcarriers,TimePositionShift+ceil(obj.PilotSpacingTime):2*obj.PilotSpacingTime:NrMCSymbols) = 1;
                  obj.PilotMatrix(FrequencyPositionShift+round(obj.PilotSpacingFrequency)    :2*obj.PilotSpacingFrequency:NrSubcarriers,TimePositionShift+1                         :2*obj.PilotSpacingTime:NrMCSymbols) = 1;
                  obj.PilotMatrix(FrequencyPositionShift+round(3/2*obj.PilotSpacingFrequency):2*obj.PilotSpacingFrequency:NrSubcarriers,(TimePositionShift+ceil(obj.PilotSpacingTime))+1:2*obj.PilotSpacingTime:NrMCSymbols) = 1;
              
              otherwise
                  error('Pilot pattern is not supported! Chose Rectangular or Diamond');
          end          
          obj.NrPilotSymbols = sum(obj.PilotMatrix(:));
          
          % preinitialize interpolation method
          switch obj.InterpolationMethod  
              case {'linear','nearest','natural'}
                  [x_pilot_pos,y_pilot_pos] = find(obj.PilotMatrix);
                  obj.InterpolationProperties = scatteredInterpolant(x_pilot_pos,y_pilot_pos,zeros(obj.NrPilotSymbols,1),obj.InterpolationMethod);         
          end              
      end 
      
      function InterpolatedChannel = ChannelInterpolation(varargin)
          obj = varargin{1};
          LSChannelEstimatesAtPilotPosition = varargin{2};
              
          switch obj.InterpolationMethod  
              case {'linear','nearest','natural'}
                  obj.InterpolationProperties.Values = LSChannelEstimatesAtPilotPosition(:);
                  [yq,xq] = meshgrid(1:size(obj.PilotMatrix,2),1:size(obj.PilotMatrix,1));
                  InterpolatedChannel = obj.InterpolationProperties(xq,yq);
              case 'FullAverage'
                  InterpolatedChannel = ones(size(obj.PilotMatrix))*mean(LSChannelEstimatesAtPilotPosition);      
              otherwise
                  error('Interpolation method not implemented');
          end
      end

      function InterpolationMatrix = GetInterpolationMatrix(varargin)
          obj = varargin{1};
          
          [x_pilot_pos,~] = find(obj.PilotMatrix);
          InterpolationMatrix = zeros(numel(obj.PilotMatrix),numel(x_pilot_pos));
          for i_pos =1:length(x_pilot_pos)
            TestDirac = zeros(size(x_pilot_pos));
            TestDirac(i_pos)=1;
            ImpulseResponse = obj.ChannelInterpolation(TestDirac);
            
            InterpolationMatrix(:,i_pos)=ImpulseResponse(:);
          end
          
      end

   end

end