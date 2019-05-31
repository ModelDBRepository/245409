%% Clear workspace
clear;

%% Add required directories for the STDP toolbox

% STDP
addpath( './Libraries/tbSTDP/' ) ;

% Sampling
addpath( './Libraries/tbSampling/' ) ;

% DoG
addpath( './Libraries/tbDoG/' ) ;

% Test-bench to explore the converged network
addpath( './Libraries/tbTestBench/' ) ;

%% Declare the project

% The project name. The ROI mapping MUST BE DONE MANUALLY later
NProject = 'Foveal';

% Create directory if it doesn't already exist
mkdir(['./Data/' NProject '/']);

% Project folder
ProjectFolder = ['./Data/' NProject '/'] ;

%% Dataset: Load paths, data about files, declare dataset PPA

% Name the dataset
Dataset.Name = 'Hibbard'; % {'KITTI','Hibbard'}

% Cortical magnification factor (CMF).
Dataset.CorticalMagnification =  1.5 ;

switch lower(Dataset.Name)
    
    case lower('KITTI')
        
        % Dataset path
        Dataset.Path = './Datasets/KITTI/' ; 
        
        % Add path to the dataset
        addpath(genpath( Dataset.Path ));
        
        % File names
        FilaNames.Folders = {'L/','R/'};
        tDir.L = dir( [Dataset.Path FileNames.Folders{1}] ) ; 
        tDir.R = dir( [Dataset.Path FileNames.Folders{2}] ) ; 
        FileNames.FileNames(:,1) = {tDir.L(~[tDir.L.isdir]).name}.' ;         
        FileNames.FileNames(:,2) = {tDir.R(~[tDir.R.isdir]).name}.' ; 
        clearvars tDir ; 
        
        % Set the number of files in the dataset
        Dataset.NFiles = size(FileNames.FileNames,1);

        % Calculating the PPA
        % Rescaling the dataset
        Dataset.Rescaling =  [1 1] ; % In pixels [horizontal x vertical]
        % Camera opening angles from the white paper on their website
        Dataset.ImageSize.Ang = [90 35]; % In degrees [horizontal x vertical]
        % Pixels on the sensor / rescaling
        Dataset.ImageSize.Pix = [1242 375]; % In pixels [horizontal x vertical]
        % PPA = Pixels / Camera angle
        Dataset.PPA = Dataset.ImageSize.Pix ./ Dataset.ImageSize.Ang ; % Approximate [horizontal x vertical]
        
    case lower('Hibbard')
        
        % Dataset path
        Dataset.Path = './Datasets/HH/' ; 
        
        % Add path to the dataset
        addpath(genpath( Dataset.Path ));
        
        % File names
        FileNames.Folders = {'L/','R/'};
        tDir.L = dir( [Dataset.Path FileNames.Folders{1}] ) ; 
        tDir.R = dir( [Dataset.Path FileNames.Folders{2}] ) ; 
        FileNames.FileNames(:,1) = {tDir.L(~[tDir.L.isdir]).name}.' ;         
        FileNames.FileNames(:,2) = {tDir.R(~[tDir.R.isdir]).name}.' ; 
        clearvars tDir ; 
        
        % Set the number of files in the dataset
        Dataset.NFiles = size(FileNames.FileNames,1);
        
        % Calculating the PPA
        % Rescaling for the dataset for faster processing
        Dataset.Rescaling =  Dataset.CorticalMagnification * [4 4] ; % In pixels [horizontal x vertical]
        % Rescaled image size in pixels = Pixels on the sensor / Rescaling
        Dataset.ImageSize.Pix = floor( [1201 1201] ./ Dataset.Rescaling ) ; % In pixels [horizontal x vertical]
        % PPA = Pixels / Camera angle after the rescaling
        Dataset.PPA = [60 60] ./  Dataset.Rescaling ; % Approximate [horizontal x vertical]
        % Camera opening angles from the white paper on their website
        Dataset.ImageSize.Ang = Dataset.ImageSize.Pix ./ Dataset.PPA ; % In degrees [horizontal x vertical]
        
end
        
%% Declare the PPA, patch-size and number of filters

PPA = Dataset.PPA; 

% Patch size: Width x Height
PatchSize.Ang =  [3 3] ; 
PatchSize.Pix = floor( PatchSize.Ang .* PPA ) ; % Width x Height

% Number of features: Number of pixels per patch x 2(ON/OFF) x 2(L/R)
NFeatures = prod(PatchSize.Pix) * 2 * 2;

%% Initialise the DoG filters and trim-size

% Define the OFF kernel size: MatlabNotation(Height x Width)
tKSizeOFF.Ang =  [1 1] ; 
tKSizeOFF.Pix = ceil( tKSizeOFF.Ang .* PPA ) ;
% Define the OFF cell centre and surround radii
tSigmaCentreSurround.OFF.Ang = [0.3 1] ;
tSigmaCentreSurround.OFF.Pix = tSigmaCentreSurround.OFF.Ang .* PPA ;

% Define the OFF cell
DoG.OFF = fCreateDoG( tKSizeOFF.Pix , PPA, ...
                      tSigmaCentreSurround.OFF.Pix , 'off');

% Spatial ON/OFF RF scaling factor
DoG.OnOffRFScalingFactor = 1 ;

% Define the ON kernel size: MatlabNotation(Height x Width)
tKSizeON.Ang = DoG.OnOffRFScalingFactor * tKSizeOFF.Ang ;
tKSizeON.Pix = ceil( tKSizeON.Ang .* PPA );
% Define the ON cell centre and surround radii
tSigmaCentreSurround.ON.Ang =  tSigmaCentreSurround.OFF.Ang ; 
tSigmaCentreSurround.ON.Pix = tSigmaCentreSurround.ON.Ang .* PPA ;

% Define the ON cell
DoG.ON = fCreateDoG( tKSizeON.Pix , PPA, ...
                     tSigmaCentreSurround.ON.Pix , 'on');

% Determine trim size (to remove the edges) after the LGN processing
TrimSize = floor ( max([DoG.ON.Size.Pix ; DoG.OFF.Size.Pix],[],1) / 2 ) ;

% Get rid of temporary variables
clearvars t* ;

%% Trimmed image size

% Image size: Width x Height
ImSize.Pix = Dataset.ImageSize.Pix - 2*TrimSize ;
ImSize.Ang = ImSize.Pix ./ PPA;

%% SAVING: DoG filtered images (beware of very large datasets)

% Option to skip
sw.SkipIt = 1 ;

if ~sw.SkipIt
    
    for nFile = 1:Dataset.NFiles
        
        % For the left and the right image
        for nLR=1:2
            
            % Read the image
            tIm = imread( [FileNames.Folders{nLR} FileNames.FileNames{nFile,nLR}] ) ;
            
            % Resize the image
            tIm = imresize(tIm, fliplr( Dataset.ImageSize.Pix ) ) ;
            
            % If the image is in colour, convert it to greyscale
            if (size(tIm,3) == 3)
                tIm = rgb2gray(imread(tIm)) ;
            end
            
            % Convert the image to double
            tIm = im2double (tIm) ;
            
            % An optional grey-world type contrast normalisation
            sw.GreyWorldNormalisation = 1; % 0 -> No, 1 -> Yes
            if sw.GreyWorldNormalisation
                % Centring each image at zero
                tIm = tIm - mean( tIm(:) ) ;
                % Scaling such that max power in any given pixel is one
                tIm = tIm / norm( tIm(:) , 2);
            end
            
            % Filtering followed by rectification. This rectification is a
            % crude method of implementing a firing threshold in the ON/OFF
            % cells - the threshold being a uniform light field.
            tFI.ON = imfilter(tIm,DoG.ON.Kernel);
            tFI.ON( tFI.ON < 0 ) = 0 ;
            tFI.OFF = imfilter(tIm,DoG.OFF.Kernel);
            tFI.OFF( tFI.OFF < 0 ) = 0 ;
            
            % Untrimmed image: MATLABFormat(height x width ) x 2(ON / OFF)
            tUntrimmedFilteredImage = cat(3, tFI.ON, tFI.OFF) ;
            
            % Trim the image to get rid of partial DoGs on the edge
            tTrimmedFilteredImage = tUntrimmedFilteredImage(...
                1+TrimSize(2) : end - TrimSize(2) , ...
                1+TrimSize(1) : end - TrimSize(1) , : );
            
            % Retinotopic RGC map: MATLABFormat(height x width ) x 2(ON /
            % OFF) x 2(L/R) x NFiles
            Retinotopic_RGCActivations(:,:,:,nLR,nFile) = tTrimmedFilteredImage ;
            
        end
    end
    
    % Save the filtered patches to HD
    save([ProjectFolder Dataset.Name '_RetinotopicRGCActivations' '.mat'],'Retinotopic_RGCActivations');
    
    
    clearvars t* n* db* Retinotopic_RGCActivations ;
    
end

%% Make ROI maps **** VERY IMPORTANT TO CHECK THIS EVERY TIME ****

% Make a position map for the patches.
[PatchMap.X, PatchMap.Y] = meshgrid(...
    linspace(-ImSize.Ang(1)/2, ImSize.Ang(1)/2, ImSize.Pix(1)),...
    linspace(-ImSize.Ang(2)/2, ImSize.Ang(2)/2, ImSize.Pix(2)));
PatchMap.Y = flipud(PatchMap.Y) ;
PatchMap.R = (PatchMap.X.^2 + PatchMap.Y.^2).^(1/2); 
PatchMap.Phi = atan2d( PatchMap.Y , PatchMap.X);

% Define edges beyond which sampling is impossible and make a mask
EdgeWidths = (ImSize.Ang / 2) - (PatchSize.Ang/2) ;
ROIMap.EdgeMask = ~( (abs(PatchMap.X) > EdgeWidths(1) ) |...
                     (abs(PatchMap.Y) > EdgeWidths(2) ) ) ;

% Manual ROI setting
ROIMap.Logical = (PatchMap.R <= 3) & ROIMap.EdgeMask ; % Foveal
             
% Cleanup
clearvars t* ;

%% Decide DoG thresholds based on a given batch

% Decide thresholds for the filters. A good strategy is to decide the
% average number of filters allowed to fire based on the threshold of the
% LGN cells.

% Load the entire dataset
load([ProjectFolder Dataset.Name '_RetinotopicRGCActivations' '.mat'],'Retinotopic_RGCActivations');

% Generate random patches per image to decide the average threshold. The
% number of samples should be as large as possible.
tNRandomSamplesPerImage = 10; % 1000
% ROI map
tParams.AvailablePosition = ROIMap.Logical ;

for nImage = 1:Dataset.NFiles
    [tRandPatch_xc, tRandPatch_yc] = ...
        fGenerateParameterSubsetPatches('Manual',tNRandomSamplesPerImage, PatchSize.Pix,ImSize, tParams ) ;
    tRF.ON (:,:,:,(1:tNRandomSamplesPerImage) + tNRandomSamplesPerImage*(nImage-1) )   =  ...
        fPickPatches( squeeze( Retinotopic_RGCActivations(:,:,1,:,nImage) ), PatchSize.Pix, tRandPatch_xc, tRandPatch_yc ) ;
    tRF.OFF(:,:,:,(1:tNRandomSamplesPerImage) + tNRandomSamplesPerImage*(nImage-1) )   =  ...
        fPickPatches( squeeze( Retinotopic_RGCActivations(:,:,2,:,nImage) ), PatchSize.Pix, tRandPatch_xc, tRandPatch_yc ) ;
end

% Decide the thresholds
FilterOutputThresholdingMethod = 'Automatic'; %{'Automatic','Graphical','Manual'}
switch lower(FilterOutputThresholdingMethod)
    
    case lower('Automatic')
        
        % The average number of ON/OFF filters allowed to fire
        DoG.AllowedToFire.ON  = floor( 0.10 * NFeatures/2) ;
        DoG.AllowedToFire.OFF = floor( 0.10 * NFeatures/2) ;
        
        % ON threshold
        [tDir.ON ,tX.ON ] = ecdf(tRF.ON(:) );
        DoG.Threshold.ON = interp1q(tDir.ON, tX.ON, 1 - DoG.AllowedToFire.ON / (NFeatures/2)) ; 
        % OFF threshold
        [tDir.OFF,tX.OFF] = ecdf(tRF.OFF(:));
        DoG.Threshold.OFF = interp1q(tDir.OFF, tX.OFF, 1 - DoG.AllowedToFire.OFF / (NFeatures/2)) ;
        
    case lower('Graphical')
        
        % Calcualte the distribution of the Gabor responses. This can be used
        % to decide the Gabor thresholds.
        
        % ON
        [tDir.ON ,tX.ON ] = ecdf(tRF.ON(:) );
        
        % OFF
        [tDir.OFF,tX.OFF] = ecdf(tRF.OFF(:));
        
        figure; hold on;
        tpDebuggerGabbies(1) = plot(tX.ON , 100* (1-tDir.ON) , 'g--' ) ;
        tpDebuggerGabbies(2) = plot(tX.OFF, 100* (1-tDir.OFF), 'r--' ) ;
        % tpDebuggerGabbies(3) = plot(tX.OnOff , (NFeatures)* (1-tF.OnOff), 'k' ) ;
        
        xlabel('Simulated threshold');
        ylabel('% of DoGs firing');
        grid minor; axis('tight');
        legend(tpDebuggerGabbies,{'ON DoGs','OFF DoGs'});
        
        disp('This option only displays the graphs. No thresholds have been set. If you want to set thresholds (perhaps based on these graphs?), use the option ''manual''.');
        
        % Set the thresholds to 0, which is the same as doing nothing to
        % the activations.
        DoG.Threshold.ON  = 0 ;
        DoG.Threshold.OFF = 0 ;    
    
    case lower('Manual')
        
        % Set the thresholds manually.
        DoG.Threshold.ON  = 1e-4 ; % Normalised 5 x 5
        DoG.Threshold.OFF = 1e-4 ; % Normalised 5 x 5

end

% Clean up temporary variables
clearvars db* t* n* h* Retinotopic_RGCActivations;

%% SAVING: Threshold and SAVE.

% Option to skip
sw.SkipIt = 1 ;

if ~sw.SkipIt
    
    % 1 <- Equal On/Off threhsolds, 0 <- Possibly unequal On/Off thresholds
    DoG.BalanceDoGOnOffThresholding = 1 ;
    if DoG.BalanceDoGOnOffThresholding
        tThreshold = min([DoG.Threshold.ON DoG.Threshold.OFF]);
        DoG.Threshold.ON = tThreshold;
        DoG.Threshold.OFF = tThreshold ;
    end
    
    % Load unthresholded retinotopic RGC map as an object. Slower but less
    % memory intensive
    objRetinotopic_RGCActivations = ...
        matfile( [ProjectFolder Dataset.Name '_RetinotopicRGCActivations' '.mat'] );
    tRetinotopic_RGCActivations = objRetinotopic_RGCActivations.Retinotopic_RGCActivations ;
    
    % Threshold file by file
    for nFile = 1:Dataset.NFiles
        
        tRF.ON = tRetinotopic_RGCActivations(:,:,1,:,nFile) ;
        tRF.ON( tRF.ON < DoG.Threshold.ON) = 0;
        Retinotopic_RGCActivations_Thresholded(:,:,1,:,nFile) = tRF.ON ;
        
        tRF.OFF = tRetinotopic_RGCActivations(:,:,2,:,nFile) ;
        tRF.OFF( tRF.OFF < DoG.Threshold.OFF) = 0;
        Retinotopic_RGCActivations_Thresholded(:,:,2,:,nFile) = tRF.OFF ;
        
    end
    
    % Save as a new variable
    save([ProjectFolder Dataset.Name '_RetinotopicRGCActivations_Thresholded' '.mat'],'Retinotopic_RGCActivations_Thresholded');
    
    % Clear the unnecessary variables
    clearvars t* n* db* obj* Retinotopic_RGCActivations_Thresholded ;
    
end

%% STDP: Declare the overall network and parameters

% Declare the overall structure
Layer = struct();

% NInputs
Layer(1).NInputs = NFeatures;
% Name
Layer(1).Name = 'L1';

% Number of neurones
Layer(1).NNeurones =  300 ; 
% Threshold
Layer(1).Threshold = 18 ;

% Initial random weights 
% {'UniformBalanced','UniformUnbalanced','NormalBalanced','NormalUnbalanced'}
WeightInitialisationScheme = 'NormalUnbalanced' ; 
switch lower(WeightInitialisationScheme)
    
    case lower('UniformUnbalanced')
        tW0_vec = rand( Layer(1).NInputs , Layer(1).NNeurones) ; 
        
    case lower('UniformBalanced')
        tW0_vec_Mono = rand( Layer(1).NInputs / 2 , Layer(1).NNeurones) ; 
        tW0_vec = [tW0_vec_Mono ; tW0_vec_Mono];    
        
    case lower('NormalUnbalanced')
        tW0_mu = 0.2 ; tW0_sd = 0.2 / 3;
        tW0_vec = tW0_sd * randn( Layer(1).NInputs , Layer(1).NNeurones ) + tW0_mu;
        tOutlierWs = ( (tW0_vec <= 0) | ( tW0_vec >= 1 ) ) ; 
        tW0_vec( tOutlierWs ) = ...
            (tW0_mu - tW0_sd) + 2*tW0_sd * rand( sum(tOutlierWs(:)) , 1 )  ; 
        
    case lower('NormalBalanced')
        tW0_mu = 0.6 ; tW0_sd = 0.4 / 2;
        tW0_vec_Mono = tW0_sd * randn(Layer(1).NInputs/2,Layer(1).NNeurones) + tW0_mu;
        tOutlierWs = ( (tW0_vec_Mono <= 0) | ( tW0_vec_Mono >= 1 ) ) ; 
        tW0_vec_Mono( tOutlierWs ) = ...
            (tW0_mu - tW0_sd) + 2*tW0_sd * rand( sum(tOutlierWs(:)) , 1 )  ; 
        tW0_vec = [tW0_vec_Mono ; tW0_vec_Mono];
        
end
Layer(1).Weight = tW0_vec;

% Inhibition strategy 
% {Uninhibited', 'LateralInhibition','LateralInhibitionCustom'}
Layer(1).InhibStrategy.Type = 'LateralInhibitionCustom'; 
% Number of neurons allowedto spike (used only with 'LateralInhibitionCustom' )
Layer(1).InhibStrategy.NFirstsNeurons = 1 ; 

% LTP and LTD rate adjustment
% Strategy used to adjust LT rates with respect to the convergence of the
% weight {'none', 'convergence'}. 'convergence' option is EXPERIMENTAL. 
Layer(1).LTRates.Strategy = 'none'; 
% Number of frames to run without updating LT rates. CAUTION: Should be
% more than the parameter LocalNeighbourhood send to fAdjustLTRates ( 5 by
% default). Anything except 'Inf' is EXPERIMENTAL.
Layer(1).LTRates.Delay = Inf;

% LTP
Layer(1).LTP.Rate = 0.005 * diag( ones(Layer(1).NNeurones,1) ) ; 
Layer(1).LTP.Bounds = 'hard'; % {'soft','hard'}
Layer(1).LTP.Mu = 0.65 ; % [0 1] Non linearity exponent
% LTD
Layer(1).LTD.Rate = - 0.75 * Layer(1).LTP.Rate ; % MUST be negative
Layer(1).LTD.Bounds = 'hard'; % {'soft','hard'}
Layer(1).LTD.Mu = 0.05; % [0 1] Non linearity exponent

% Clean up temporary variables
clearvars db* t* n* h* ;

%% STDP: Run L1 to convergence, with lateral inhibition

% Number of patches to test
Training.NPatches = 1e5 ;
% Number of batches
Training.NMaxPatchesPerBatch = 100 ; 
% ROI Map. Temporary variable to store parameters.
tParams.AvailablePosition = ROIMap.Logical ;
% Counter for iterations
nTotalPatchesRun = 0;

% Load mat file of thresholded RGC outputs
tRGC_Object = matfile([ProjectFolder Dataset.Name '_RetinotopicRGCActivations_Thresholded' '.mat']);
tRGC_Activations = tRGC_Object.Retinotopic_RGCActivations_Thresholded ;

% Total processing time
tTotalProcessingTime = 0;

% Time lapse weights file (to store the time-course of weights)
tFile_WsTimeLapse = ['Ws_TimeLapse_' Dataset.Name '_' datestr(now,'ddmmyyHHMMSS') '.mat'] ;

% Parameters for snapshots of weights
tNSnapshotsForWsTimeLapse = 200 + 1 ; % 1 for the random initial weights
tApproxPatchesBetweenSnapshot = round( Training.NPatches / tNSnapshotsForWsTimeLapse ) ;

% Initialise the .mat file for weight snapshots
dbSnapshotWeights = 1;
if dbSnapshotWeights
   
    % Counter for patches run between snapshots
    nPatchesSinceLastSnapshot = 0;
    
    % Declare a dummy first frame so that matlab doesn't get confused
    Ws_TimeLapse(:,:,1) = zeros( size(Layer(1).Weight) ) ;
    
    % Store the random initial weights
    Ws_TimeLapse(:,:,2) = Layer(1).Weight ;
    % fplotUpdateWeight(Ws_TimeLapse(:,:,1) ,1) ;
    
    % Save weights
    save([ProjectFolder tFile_WsTimeLapse],'Ws_TimeLapse','-v7.3');
    
    % clear the workspace variable
    clearvars Ws_TimeLapse ;
end

% Counter for batches
nBatchN = 0;

% Keep running the network till Trainin.NPatches patches are run
while (nTotalPatchesRun < Training.NPatches)
    
    % Increment batch number
    nBatchN = nBatchN + 1;
    
    % Number of patches to run in this batch
    if ( (Training.NPatches - nTotalPatchesRun) >= Training.NMaxPatchesPerBatch )
        nPatchesToRunThisBatch = Training.NMaxPatchesPerBatch ;
    else
        nPatchesToRunThisBatch = (Training.NPatches - nTotalPatchesRun);
    end
    
    % Total number of patches run
    nTotalPatchesRun = nTotalPatchesRun + nPatchesToRunThisBatch ;
        
    % Initialise random images to sample from
    nImagesToSample = randi([1 Dataset.NFiles],nPatchesToRunThisBatch,1) ;
    
    % Initialise random sampling positions within the ROI
    [tRandPatch_xc, tRandPatch_yc] = ...
        fGenerateParameterSubsetPatches('Manual',nPatchesToRunThisBatch, PatchSize.Pix, ImSize, tParams ) ;
    
    % Randomly sample patches for the current batch
    for nPatchBeingSampled = 1:nPatchesToRunThisBatch
        
        tFilterActivations(:,:,:,:,nPatchBeingSampled) = ...
            fPickPatches( squeeze(tRGC_Activations(:,:,:,:,nImagesToSample(nPatchBeingSampled)) ) , ...
                          PatchSize.Pix, tRandPatch_xc(nPatchBeingSampled), tRandPatch_yc(nPatchBeingSampled) ) ;
        
    end
    
    % Convert the filter intensity vector to a latency vector
    tFilterLatencies_vec = reshape( fIntencity2Latency(tFilterActivations,'1/x') , [] , nPatchesToRunThisBatch) ;

    % Clear the batch-activations
    clearvars tFilterActivations ;
    
    % Run the network    
    Layer = fmultiLayersSTDP(Layer , tFilterLatencies_vec, 1, 0, 0);
    
    % Plot the last batch for debugging
    dbPlotWeightsLastBatch = 0 ;
    if ( (dbPlotWeightsLastBatch) && ( nTotalPatchesRun == Training.NPatches ) )
        % Only plot the first 20 neurones
        fplotUpdateWeight ( Layer(1).Weight(:,1:20) , 1) ;
    end    
    
    % Update processing time
    tTotalProcessingTime = tTotalProcessingTime + Layer(1).Diagnostic.ProcessingTime ;
    
    % Display processing times per batch for diagnostics
    dbDisplayProcessingTimesPerBatch = 1 ;
    if dbDisplayProcessingTimesPerBatch
        disp(['Convergence ' Layer(1).Name ': Batch# ' num2str(nBatchN) ', ' ...
            num2str(nPatchesToRunThisBatch) ' patches,  ' ...
            num2str(Layer(1).Diagnostic.ProcessingTime) ' s.']);
    end
    
    % Store weights after each run for debugging
    if dbSnapshotWeights 
        % Increment patches since last snapshot
        nPatchesSinceLastSnapshot = nPatchesSinceLastSnapshot + nPatchesToRunThisBatch ;
        
        % If enough patches processed, store weights-snapshot
        if (nPatchesSinceLastSnapshot >= tApproxPatchesBetweenSnapshot) 
            
            % Load the time lapse file and store the weights
            tWs_TimeLapse = matfile([ProjectFolder tFile_WsTimeLapse],'Writable',true);
            %         tWs_TimeLapse.Ws_TimeLapse(:,:,nBatch + 1) = Layer(1).Weight ;  
            tWs_TimeLapse.Ws_TimeLapse(:,:,end + 1) = Layer(1).Weight ;  
            
            % Reset the patch counter
            nPatchesSinceLastSnapshot = 0;

        end
    end
    
end

% Display overall processing times
disp(['Convergence ' Layer(1).Name ': Total ' num2str(nBatchN) ' batch(es), ' ...
    num2str( Training.NPatches ) ' patches,  ' num2str(tTotalProcessingTime) ' s.']);

% Clean up temporary variables
clearvars db* t* n* h* ;

%% Calculate converged L1 receptive fields

% Algorithm for RF calculation
Layer(1).RFs.Algorithm = 'LinearApproximation';
% Backpropagated weights
Layer(1).RFs.WsBP_vec = Layer(1).Weight ;

% RFs
Layer(1).RFs.RFs = fApproximateDoGRF( Layer(1).RFs.WsBP_vec , ...
                                      PatchSize.Pix   , ...
                                      DoG             , ...
                                      Layer(1).RFs.Algorithm );

%% Fit gabors to the converged neurones

% Option to skip
sw.SkipIt = 0 ;

if ~sw.SkipIt    
    
% Sinusoid and Gaussian orientation constraint
sw.SinusoidGaussianOrientationEqual = 1 ;

% Define temporary variables to make the parfor loops faster
tNeuronesToFit = 1:Layer(1).NNeurones ;
tRFs = Layer(1).RFs.RFs ;
tLayerName = Layer(1).Name ;
tSinusoidGaussianOrientationEqual = sw.SinusoidGaussianOrientationEqual ; 

% Start the timer
cTotalFitting = tic; 

% Use the parallel for-ce, Luke
parfor nnNeurone = tNeuronesToFit 
    cSingleFit = tic;
    
    nNeurone = tNeuronesToFit(nnNeurone) ;
    
        tFittedGabors(nnNeurone) = fFitGabors( tRFs(:,:,:,nNeurone) , PPA , 0,...
        sw.SinusoidGaussianOrientationEqual, 3 ) ; 
    
    ttProcessingTime = toc(cSingleFit);

    % Display processing times per batch for diagnostics
    dbDisplayProcessingTimesPerBatch = 1;
    if dbDisplayProcessingTimesPerBatch
        disp(['Gabor fitting ' tLayerName ': Neurone# ' num2str(nNeurone) ...
              ', ' num2str(ttProcessingTime) ' s.']);
    end
    
end

% Total processing time for the fitting process
tTotalProcessingTime = toc(cTotalFitting) ; 

% Gather output of the parfor
Layer(1).FittedGabors = tNeuronesToFit ;
for nnNeurone = 1:numel(tNeuronesToFit) 
    nNeurone = tNeuronesToFit(nnNeurone) ;
    
    Layer(1).FittedGaborParameters_vecs(1:9,:,nNeurone) = tFittedGabors(nnNeurone).Parameters_vec ;
    Layer(1).FittedGaborParameters_vecs(10 ,:,nNeurone) = tFittedGabors(nnNeurone).R2.' ;
    
end


% Display overall processing times
disp(['Gabor fitting ' Layer(1).Name ': All neurone(s) ' num2str(numel(tNeuronesToFit)) ...
    ', ' num2str(tTotalProcessingTime) ' s.']);

% Clean up temporary variables
clearvars db* t* n* h* ;

end

%% SAVING: The workspace up to this point

% Option to skip
sw.SkipIt = 0 ;

if ~sw.SkipIt 
    
    save([ProjectFolder 'ws_' Dataset.Name '_' datestr(now,'ddmmyyHHMMSS') '.mat']);
    
end

%% Fit Horizontal DTCs using binocular correlation method

% Option to skip
sw.SkipIt = 0 ;

if ~sw.SkipIt
    
    % Create a field for Binocular correlation DTCs
    Layer(1).DTCs_BC.NeuronesTested = 1:Layer(1).NNeurones ;
    
    % Calculate the BC DTCs
    for nnNeurone = 1:numel(Layer(1).DTCs_BC.NeuronesTested)
        nNeurone = Layer(1).DTCs_BC.NeuronesTested( nnNeurone );
        
        % Extract the RF
        tRF = Layer(1).RFs.RFs(:,:,:,nNeurone) ;
        % Calculat the horizontal binocular correlation
        tDTC = fHorizontalDTC_Correlation( tRF, PPA, PatchSize ) ;
        
        Layer(1).DTCs_BC.DTCs_raw(:,nnNeurone) = tDTC.DTC_raw ;
        Layer(1).DTCs_BC.GaborFit_vec(:,nnNeurone) = tDTC.GaborFit_vec ;
        Layer(1).DTCs_BC.DTCs_Gabor(:,nnNeurone) = tDTC.GaborFit ;
        if nnNeurone == 1
            Layer(1).DTCs_BC.Disparities = tDTC.DisparityRange.' ;
        end
    end
    
    % Save the results with BC DTCs
    sw.SaveIt = 1;
    if sw.SaveIt
        
        save([ProjectFolder 'ws_DTCs_BC_' Dataset.Name '_' datestr(now,'ddmmyyHHMMSS') '.mat']);
        
    end
            
end
