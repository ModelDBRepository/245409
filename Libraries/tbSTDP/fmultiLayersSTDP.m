function Layers = fmultiLayersSTDP(Layers, InitialInputs, UdpdateWeight, SaveFiringLatencies, Plot, SaveDiagnostic)
% The purpose of this fuction is to simplify use of network with multiple layers. It is still in an EXPERIMENTAL stage, but seems to be quite stable. 
% 
% For example:      
% To compute output and update weight for layer 1, use firing latencies of layer 1 to process layer 2, etc.. until the last layer, :
% Layers = fmultiLayersSTDP(Layers, InitialInputs, [1 1 (..) 1], [1 1 (..) 1]);
%
% % % % % % % % % % % % % % % % % % % % % % % % % % 
%       INPUTs
% % % % % % % % % % % % % % % % % % % % % % % % % % 
%       Layers(i).                          <struct>
% Layers(i).Weight                          <matrix> Synapse weight, specific to each neurons and each inputs [Nn*Ni]
% Layers(i).NNeurones                       [1] The number of neurones
% Layers(i).Threshold                       [1] Threshold for the current layer
% Layers(i).Name                            <string> Layer names ('L1', 'L2', ...)
% Layers(i).InhibStrategy.Type              <string> Strategy of ihibition : 'Uninhibited', 'LateralInhibition',  'LateralInhibitionCustom'
% Layers(i).InhibStrategy.NFirstsNeurons    [1] Numbers of neurons allowedto spike (used only with 'LateralInhibitionCustom' )
% Layer(i).LTRates.LTP                      <matrix> Determine initial rate of the "Long Term Potentialization"
% Layers(i).LTRates.LTD                     <matrix> Determine initial rate of the "Long Term Depression"
% Layers(i).LTRates.Strategy                <string> Strategy used to adjust LT rates w.r.t. the convergence of the weights (Available: 'none', 'convergence'). This feature is EXTREMELY EXPERIMENTAL at this stage.
% Layers(i).LTRates.Delay                   [1] Number of frames without updates LTrates, CAUTION : Should be more than the parameter LocalNeighbourhood sent to fAdjustLTRates ( 5 by default)
%      
% InitialInputs                             <matrix> Spike latencies from each elementary input, and for each frame [Ni*Nframe] (Rq: only latencies here, no structure)
% 
% UdpdateWeight                             <boolean vector> Determine if weight of layer 'i', need to be updated ( Should layer i lear or not ? )
% SaveFiringLatencies                       <boolean vector> Determine if output of layer 'i', need to be saved (CAUTION, if not, the next layer can't be computed..)
% Plot                                      <boolean vector> Determine if fSTDP should call fplotUpdateWeight (false by default)
% SaveDiagnostic                            <boolean vector> Determine if Diagnostic.WeightOverTime & Diagnostic.WeightOverTime  need to be saved (false by default)
% INTEGRER option LayerbyLayer ... 1st Layer until convergence, then second Layer, etc... 
% % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % 
%       OUTPUTs
% 
%       Layers(i).                          <struct> 
% Layers(i).Weight                          <matrix> Contain the new weight for layer with active weight's udpdate   
% Layers(i).FiringLatencies                 <matrix> Firing latencies of the layer (only for layer with active output)
% Layers(i).Diagnostic                      <struct> 
% ( Layers(i). ... )                        Other subfields listen above are not changed
% 
% RQ Changer .LTRates. ???


%% Internal parameter. Layers expected to be binocular
binocularLayer = 1 ; 

%% Check INPUTS
NLayer = size(Layers,2); 
if NLayer ~= size(UdpdateWeight,2) || NLayer ~= size( SaveFiringLatencies,2)
    disp('ERROR from fmultiplayerSTDP: The number of layer doesn''t match the size of "SaveFiringLatencies" and/or "UdpdateWeight" ');
    
else

%% Default parameters
if nargin < 6
   SaveDiagnostic = zeros( 1, size(UdpdateWeight,2) ); 
end
if nargin < 5
   Plot = zeros( 1, size(UdpdateWeight,2) ); 
end

%% Initialise temp variable nextInput
nextInput.Latencies  = InitialInputs; 
nextInput.Nframe = size( InitialInputs, 2 );

for i = 1:NLayer

    if (UdpdateWeight(i)) == 1 && (SaveFiringLatencies(i) == 1)
        % Set param 
        Param.UpdateWeightOnOff = 1;
        Param.SaveLatencies = 1;
        if Plot(i)
             Param.Plot.OnOff = 1;  
        else
            Param.Plot.OnOff = 0;
        end     
        
        if sum(i == binocularLayer)
            Param.Plot.Binocular =1;
        else
            Param.Plot.Binocular =0;
        end
          
        if SaveDiagnostic
            Param.SaveDiagnostic = true;
        else
            Param.SaveDiagnostic = false;
        end
        
        [ Layers(i).Weight, Layers(i).FiringLatencies, Layers(i).Diagnostic ] = fSTDP (Layers(i), nextInput, Param);

        
        nextInput.Latencies = Layers(i).FiringLatencies ; 
        nextInput.Nframe = size(Layers(i).FiringLatencies, 2 );
        
    elseif  UdpdateWeight(i) == 1 && SaveFiringLatencies(i) ==0
         % Set param
        Param.UpdateWeightOnOff = 1;
        Param.SaveLatencies = 0;
        if Plot(i)
            Param.Plot.OnOff = 1;
        else
            Param.Plot.OnOff = 0;
        end      
        
        if sum(i == binocularLayer)
            Param.Plot.Binocular =1;
        else
            Param.Plot.Binocular =0;
        end      
        
        if SaveDiagnostic
            Param.SaveDiagnostic = true;
        else
            Param.SaveDiagnostic = false;
        end
        
        [ Layers(i).Weight, ~ , Layers(i).Diagnostic ] = fSTDP (Layers(i), nextInput, Param);
    
        nextInput.Latencies = NaN;
        nextInput.Nframe = NaN;
        
    elseif  UdpdateWeight(i) == 0 && SaveFiringLatencies(i) ==1        
         % Set param
        Param.UpdateWeightOnOff = 0;
        Param.SaveLatencies = 1;
        
        if Plot(i)
            Param.Plot.OnOff = 1;
        else
            Param.Plot.OnOff = 0;
        end   
        if sum(i == binocularLayer) 
            Param.Plot.Binocular =1;
        else
            Param.Plot.Binocular =0;
        end          
        
        if SaveDiagnostic
            Param.SaveDiagnostic = true;
        else
            Param.SaveDiagnostic = false;
        end
        
        [ ~, Layers(i).FiringLatencies, Layers(i).Diagnostic ] = fSTDP (Layers(i), nextInput, Param);
        
        nextInput.Latencies = (Layers(i).FiringLatencies) ; 
        nextInput.Nframe = size(Layers(i).FiringLatencies, 2 );
        
    end
    
end

end

% Cleanup
clearvars -except Layers ;
