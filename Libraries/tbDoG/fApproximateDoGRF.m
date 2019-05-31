function RFs = fApproximateDoGRF( Weights_vec, PatchSizePix, DoGKernels, RFDisplayAlgorithm )
% Given a set of weights and bases, this function returns the RF.
% 
% Inputs:
% Weights               : [NFilters NNeurones] Weights. 
% PatchSizePix          : [Width Height] Receptive field size in pixels
% DoGKernels            : Structures with the DoG information 
% RFDisplayAlgorithm	: {'Weights','WeightConvolutions','LinearApproximation'} Algorithm to use
%                         for calculating the RFs
%
% Outputs: 
% RFs                   : [KSize NLR NNeurones]

%% Extract the relevant parameters

NNeurones = size( Weights_vec , 2) ; 
DoG = DoGKernels ;

%% Reshape the weights using the patch size

% Reshape the weights of the layer 1 neurones. To understsand the
% dimension, follow the size of FI ! 
Ws = reshape( Weights_vec , ...
              PatchSizePix(2), PatchSizePix(1), ...
              2, 2, ... 
              NNeurones );

%% Calculate RFs

RFs = zeros( PatchSizePix(2) , PatchSizePix(1) , 2 , NNeurones );

% Calcualte RFs based on the algorithm
switch lower(RFDisplayAlgorithm)
    
    case lower('Weights')
        % Simple WOn - WOff
        RFs = squeeze( Ws(:,:,1,:,:) - Ws(:,:,2,:,:) ) ; 
        
    case lower('WeightConvolutions')

        RFs  = imfilter( squeeze( Ws(:,:,1,:,:) ),DoG.ON.Kernel)  + ... 
               imfilter( squeeze( Ws(:,:,2,:,:) ),DoG.OFF.Kernel) ;
    
    case lower('LinearApproximation')
        
        tTrims.ON  = floor( DoG.ON.Size.Pix / 2 );
        tTrims.OFF = floor( DoG.OFF.Size.Pix / 2 );
        
        % Declare a blank RF with the trim included
        tBlankRF.ON  = zeros( PatchSizePix(2) + 2*tTrims.ON(2), ...
                              PatchSizePix(1) + 2*tTrims.ON(1) ) ;
        tBlankRF.OFF = zeros( PatchSizePix(2) + 2*tTrims.OFF(2), ...
                              PatchSizePix(1) + 2*tTrims.OFF(1) ) ;
                         
        % Reconstruct the RF linearly.
        for nDx = 1:PatchSizePix(1)
            for nDy = 1:PatchSizePix(2)
                    for nLR = 1:2 
                        for nNeurone = 1:NNeurones
                            
                            tRF.ON = tBlankRF.ON;
                            tRF.ON( nDy -1 + ( 1:DoG.ON.Size.Pix(2) ) , ...
                                    nDx -1 + ( 1:DoG.ON.Size.Pix(1) ) ) = ...
                                DoG.ON.Kernel * Ws(nDy,nDx,1,nLR,nNeurone) ;
                            
                            tRF.OFF = tBlankRF.OFF;
                            tRF.OFF( nDy -1 + ( 1:DoG.OFF.Size.Pix(2) ) , ...
                                     nDx -1 + ( 1:DoG.OFF.Size.Pix(1) ) ) = ...
                                DoG.OFF.Kernel * Ws(nDy,nDx,2,nLR,nNeurone) ;
                            
                            RFs(:,:,nLR,nNeurone) = ...
                                RFs(:,:,nLR,nNeurone) + ...
                                tRF.ON ( tTrims.ON(2)  + ( 1:PatchSizePix(2) ) ,   ...
                                         tTrims.ON(1)  + ( 1:PatchSizePix(1) ) ) + ...
                                tRF.OFF( tTrims.OFF(2) + ( 1:PatchSizePix(2) ) ,   ...
                                         tTrims.OFF(1) + ( 1:PatchSizePix(1) ) ) ;
                        end
                    end
            end
        end
        
end

