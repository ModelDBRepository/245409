function [Weight, FiringLatencies, Diagnostic] = fSTDP (Layer, TrainingSet, Param)
% This function is our STDP implementation.
%
% % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % %
%       INPUTs
%%%
% TrainingSet.Latencies                 [Ni*Nframe] Spike latencies from each elementary input, and for each frame
% TrainingSet.Ni                        [1] The number of elementary inputs, current impletation assumes each input connects to all neurones.
% TrainingSet.Nframe                    [1] The frame-number, our abstract time
%
% Layer.Weight                          [NNeurones x Ni] Synapse weight, specific to each neurone and each input
% Layer.NNeurones                       [1] The number of neurones
% Layer.Threshold                       [1] Threshold for the current layer
% Layer.Name                            <string> Layer name ('L1', 'L2', ..)
% Layer.InhibStrategy.Type              <string> Strategy of ihibition : 'Uninhibited', 'LateralInhibition',  'LateralInhibitionCustom'
% Layer.InhibStrategy.NFirstsNeurons    [1] Numbers of neurons allowed to spike (used only with 'LateralInhibitionCustom' ) (1/100 of Neurons by default)
% Layer.LTP.Bounds                      <string> {'soft','hard'} Bounding method. A soft bound uses x*(1-x) as a bound for the potentiation whereas a hard bound strictly implements a (0,1) range for the final weights.
% Layer.LTP.Rate                        [NNeurones x NNeurones] Determine initial LTP rate
% Layer.LTP.Mu                          [1] Non linearity exponent.
% Layer.LTD.Bounds                      <string> {'soft','hard'} Bounding method. A soft bound uses x*(1-x) as a bound for depression whereas a hard bound strictly implements a (0,1) range for the final weights.
% Layer.LTD.Rate                        [NNeurones x NNeurones] Determine initial LTDP rate
% Layer.LTD.Mu                          [1] Non linearity exponent.
% Layer.LTRates.Strategy                <string> Strategy used to adjust LT rates (Available : 'none', 'convergence')
% Layer.LTRates.Delay                   [1] Number of frames without updates LTrates, CAUTION : should be more than parameter LocalNeighbourhood send to fAdjustLTRates ( 5 by default)
%
% Param.Plot.OnOff                      <boolean> Plot schematic view of the layer ?
% Param.Plot.Binocular                  <boolean> Is the layer binocular ?
% Param.UpdateWeightOnOff               <boolean> Update synaptic weights ?
% Param.SaveDiagnostic                  <boolean> Save Diagnostic structure "WeightOverTime" & "ConvergenceIndice" ? (false by default to save RAM)
% Param.SaveLatencies                   <boolean>
% % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % %
%       OUTPUTs
%
% Weight                        [Nsynapses x NNeurones] Final weights
% FiringLatencies               [NNeurones x Nframes] Firing Latencies, INF for non spiking neurons
% Diagnostic.Error              <string> Contains error message
% Diagnostic.ConvergenceIndice  [Nframes x NNeurones]
% Diagnostic.WeightOverTime     [Nsynapses x NNeurones x Nframes]
% Diagnostic.TimeProcess        [1]
%

%% Start the clock
TStart__ = tic;

%% Check Inputs
Diagnostic.Error = 'No error reported';
% Weights should be positives
if logical(sum(Layer.Weight(:) < 0))
    disp(['WARNING from fSTDP: the layer ' Layer.Name ' has negative weight(s) ! ']);
    Diagnostic.Error = 'Negative Weight !';
end
% Input should match with Weight
if size(Layer.Weight,1) ~= size(TrainingSet.Latencies,1)
    disp(['ERROR from fSTDP: the number of inputs given to layer ' Layer.Name ' doesnt match the number of synapses ! (Weight) ']);
    disp(' size(Layer.Weight,2) ~= size(TrainingSet.Latencies,1)');
    Diagnostic.Error = 'Input and Weight have different sizes !';
end

% Input should match with Nframe
if TrainingSet.Nframe ~= size(TrainingSet.Latencies,2)
    disp(['ERROR from fSTDP: the number of inputs given to layer ' Layer.Name ' doesnt match the number of Nframes ! ']);
    disp(' TrainingSet.Nframe ~= size(TrainingSet.Latencies,2) ');
    Diagnostic.Error = 'Input and Nframe have different sizes !';
end

if (~Param.SaveDiagnostic)  && Param.Plot.OnOff
    disp('ERROR FROM fSTDP: Currently, you can not Plot evolution of the weight with Param.SaveDiagnostic OFF ...')
    disp('Param.SaveDiagnostic is False  && Param.Plot.OnOff is ON, please deactive Plot option or active SaveDiagnostic.. (or wait for futur update)');
end

%% Defaults values

% Assign default value to NFirstsNeurons :
if isfield(Layer.InhibStrategy, 'NFirstsNeurons')
    NFirstsNeurons =  Layer.InhibStrategy.NFirstsNeurons;
else
    NFirstsNeurons = ceil(Layer.NNeurones/100);
end

% Default delay used before Adjust LTRates
if ~isfield( Layer.LTRates,'Delay' )
    Layer.LTRates.Delay = 5;
end

if ~isfield( Param,'SaveDiagnostic' )
    Param.SaveDiagnostic  = false;  
end

%% STDP implementation

NNeurones = Layer.NNeurones;

%% Initialize outputs's variables
if Param.SaveDiagnostic
    ConvergenceIndices = zeros(TrainingSet.Nframe,NNeurones);
    Diagnostic.WeightOverTime = zeros([size(Layer.Weight ) TrainingSet.Nframe]);
end

if Param.SaveLatencies
    FiringLatencies  = Inf*ones(TrainingSet.Nframe,NNeurones) ;
    % Heuristics for tuning the network. Infs code no activity. Only stored for
    % the batch that is being processed.
else
    FiringLatencies  = NaN;
end

TotalNumberOfSpikingNeurons = 0;
NumberOfIterationWithNoSpikingNeurons = 0;

%% Process the input data batch
for indT = 1:TrainingSet.Nframe
    
    % Serial input from the patches
    tPreSynapticLatencies = TrainingSet.Latencies(:,indT);
    
    tPresynapticFireNoFire = ~isinf(tPreSynapticLatencies);
    
    [~,indSortedPresynapticLatencies] = sort(  tPreSynapticLatencies ) ;
    
    projectedMembranePotentials = ...
        cumsum( ...
        ( tPresynapticFireNoFire(indSortedPresynapticLatencies,:) * ones(1,NNeurones) ) .*  ...
        Layer.Weight (indSortedPresynapticLatencies,:) ,...
        1 );
    
    indSpikingNeurones = find( projectedMembranePotentials(end,:) >= Layer.Threshold ) ;
    
    indLatencyMaxLTP = 1 + sum( (projectedMembranePotentials(:,indSpikingNeurones) < Layer.Threshold) ,1 ) ;
    
    if isempty(indSpikingNeurones)
        NumberOfIterationWithNoSpikingNeurons = NumberOfIterationWithNoSpikingNeurons + 1;
    end
    
    if ( ~isempty(indSpikingNeurones) && Param.UpdateWeightOnOff  )
        
        if (~strcmp(Layer.LTRates.Strategy, 'none')) && indT > Layer.LTRates.Delay    % for the first x inputs do nothing
            
            [Layer.LTP, Layer.LTD]  = fAdjustLTRates(Layer.LTP, Layer.LTD, Layer.LTRates, ConvergenceIndices, indT);
            
        end
        
        switch lower(Layer.InhibStrategy.Type )
            
            case lower('Uninhibited')
                % No lateral interaction
                
                Masks.LTP = zeros(size(Layer.Weight ));
                for nSpikingNeurone = 1:numel(indSpikingNeurones)
                    Masks.LTP(...
                        indSortedPresynapticLatencies( 1:indLatencyMaxLTP(nSpikingNeurone) ), ...
                        indSpikingNeurones(nSpikingNeurone) ...
                        ) = 1;
                end
                
                Masks.LTD = zeros(size(Layer.Weight ));
                Masks.LTD(:,indSpikingNeurones ) = ~ Masks.LTP(:,indSpikingNeurones ) ;
                
            case lower('LateralInhibition')
                % Lateral inhibition
                
                [~,indNeuroneFiringRanks] = ...
                    sort( indLatencyMaxLTP + rand( size(indLatencyMaxLTP) ) - 0.5,...
                    'ascend') ;
                indFirstNeurone = indNeuroneFiringRanks(1) ;
                
                indLatencyMaxLTP = indLatencyMaxLTP(indFirstNeurone);
                indSpikingNeurones = indSpikingNeurones(indFirstNeurone);
                
                Masks.LTP = zeros(size(Layer.Weight ));
                Masks.LTP(...
                    indSortedPresynapticLatencies( 1:indLatencyMaxLTP ), ...
                    indSpikingNeurones ...
                    ) = 1;

                Masks.LTD = Masks.LTP;
                Masks.LTD(:,indSpikingNeurones ) = ...
                    ~Masks.LTD(:,indSpikingNeurones ) ;
                
            case lower('LateralInhibitionCustom')
                % Lateral inhibition Custom,
                
                NFirstsNeuronstemp = NFirstsNeurons;
                
               
                if  sum(logical(indSpikingNeurones)) < NFirstsNeurons
                    NFirstsNeuronstemp = sum(logical(indSpikingNeurones));
                end
                
                [~,indNeuroneFiringRanks] = ...
                    sort( indLatencyMaxLTP + rand( size(indLatencyMaxLTP) ) - 0.5,...
                    'ascend') ;
                
                % Only consider the spikes for the 'NFirstsNeurons' which spikes
                indFirstNeurones = indNeuroneFiringRanks(1:NFirstsNeuronstemp) ;
                indSpikingNeurones = indSpikingNeurones(indFirstNeurones);
                indLatencyMaxLTP = indLatencyMaxLTP(indFirstNeurones);
                
                Masks.LTP = zeros(size(Layer.Weight ));
                for nSpikingNeurone = 1:numel(indSpikingNeurones)
                    Masks.LTP(...
                        indSortedPresynapticLatencies( 1:indLatencyMaxLTP(nSpikingNeurone) ), ...
                        indSpikingNeurones(nSpikingNeurone) ...
                        ) = 1;
                end
                
                Masks.LTD = Masks.LTP;
                Masks.LTD(:,indSpikingNeurones ) = ...
                    ~Masks.LTD(:,indSpikingNeurones ) ;
                
                TotalNumberOfSpikingNeurons = TotalNumberOfSpikingNeurons + NFirstsNeuronstemp  ;
            otherwise
                disp('ERROR !! from function fSTDP, the value of the argument "Layer.InhibStrategy.Type" is incorrect ');
                Diagnostic.Error = 'Invalid argument in "Layer.InhibStrategy.Type" ';
                
        end
        
        %% update Weight
        
        Layer.Weight = fUpdateWeights(Layer.Weight, Masks, Layer.LTP, Layer.LTD );
        
    end
    
    % Update the latencies matrix
    if Param.SaveLatencies
        FiringLatencies(indT,indSpikingNeurones)  = indLatencyMaxLTP;
    end
    
    if Param.SaveDiagnostic
        ConvergenceIndices(indT,:) = mean( abs( Layer.Weight  - (Layer.Weight  > 0.5) ) , 1 )  ;
        
        % Save the weights
        Diagnostic.WeightOverTime(:,:,indT) = Layer.Weight ;
    end
    
end

switch lower(Layer.InhibStrategy.Type )
    case lower('LateralInhibitionCustom')
        disp( ['NO SPIKE ITR: ' ,num2str(NumberOfIterationWithNoSpikingNeurons)   ,'; AVG(SPIKING ONLY): ', num2str(TotalNumberOfSpikingNeurons/(TrainingSet.Nframe - NumberOfIterationWithNoSpikingNeurons )), '; MAX SPIKE: ',num2str(NFirstsNeurons)] ) ;
    otherwise
        disp( ['NO SPIKE ITR: ' ,num2str(NumberOfIterationWithNoSpikingNeurons) ]);
        
end

Weight = Layer.Weight ;

if Param.SaveLatencies
    FiringLatencies = FiringLatencies.';
end
if Param.SaveDiagnostic
    Diagnostic.ConvergenceIndice = ConvergenceIndices;
end
%% Plot

if Param.Plot.OnOff
    
    fplotUpdateWeight( Diagnostic.WeightOverTime, Param.Plot.Binocular );
    
end


%% Clean outputs..

Diagnostic.ProcessingTime = toc(TStart__);
clearvars -except Weight FiringLatencies Diagnostic ;
