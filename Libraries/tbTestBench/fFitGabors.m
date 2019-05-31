function FittedGabor = fFitGabors( RFs, PPA, PlotFittedGabor01, ForceOrientationDependance, BestOfNRuns )
% Function for fitting Gabors to RFs. 
% 
% Inputs:
% RFs: Receptive fields
% PPA: Pixels per angle
% PlotFittedGabor01: Plot the fit or not
% ForceOrientationDependance: Make Gabor envelope and sinusoid have the same orientation or not
% BestOfNRuns: Take the best fit out of N runs of the fitting function
% 
% Outputs: 
% FittedGabor: The fitted Gabor!

%% Extract parameters from the input

% Determine if it is a binocular or a monocular fit
NLR = size(RFs,3);

% Determine the size of the grid
PatchSize.Pix = [ size(RFs,2) size(RFs,1) ] ;
PatchSize.Ang = PatchSize.Pix ./ PPA ;

% Create a mesh
[Mesh.XX, Mesh.YY] = meshgrid(...
    linspace(-PatchSize.Ang(1)/2, PatchSize.Ang(1)/2, PatchSize.Pix(1)),...
    linspace(-PatchSize.Ang(2)/2, PatchSize.Ang(2)/2, PatchSize.Pix(2)));
Mesh.YY = flipud(Mesh.YY) ;

% Extract a measure of the total variance from the data
for nLR = 1:NLR
    % Variance estimate (unbiased estimator)
    SS_tot(nLR) = var( reshape( RFs(:,:,nLR),[],1) , 0 , 1 ) ; 
end

%% Decide the starting point for the optimisation

% A meshgrid of frequencies
[F.XX, F.YY] = ...
    meshgrid( -0.5*PPA(1) : PPA(1)/PatchSize.Pix(1) : 0.5*PPA(1) - PPA(1)/PatchSize.Pix(1) , ...
              -0.5*PPA(2) : PPA(2)/PatchSize.Pix(2) : 0.5*PPA(2) - PPA(2)/PatchSize.Pix(2) );

for nLR = 1:2
    
    % Monocular field
    RF = RFs(:,:,nLR);
    
    % Find max indices for the spatio-temporal domain. 
    [RFMax, tInd_RFMax] = max( reshape( RF ,1,[]) );
    [~, tInd_RFMin] = min( reshape( RF ,1,[]) );    
    
    % Gabor scalings
    Gabor.ScalingFactor(nLR) = RFMax ;
    
    % The frequency and orientation.
    RFf.FFT = fftshift( fft2( flipud( RF ) ) ) ;
    RFf.Abs = abs( RFf.FFT); 
    RFf.Phase = rad2deg( angle( RFf.FFT) );
    [~, tInd_FMax] = max( reshape( RFf.Abs ,1,[]) );
    
    
    % Gabor Sinusoid frequency is f_Gabor = max_f{ |FFT(f)| }
    Gabor.Sinusoid.Frequency(nLR) = norm([F.XX(tInd_FMax) F.YY(tInd_FMax)]); 
    % Gabor sinusoid orientation is Theta_Gabor = arg( f_opt ) in [0,180)
    Gabor.Sinusoid.Theta(nLR) = atan2d( F.YY(tInd_FMax), F.XX(tInd_FMax) );
    if ( Gabor.Sinusoid.Theta(nLR) < 0 )
        Gabor.Sinusoid.Theta(nLR) = Gabor.Sinusoid.Theta(nLR) + 180 ;
    end
    
    % Gabor centre
    Gabor.Centre(:,nLR) = [mean(Mesh.XX([tInd_RFMin tInd_RFMax])) mean(Mesh.YY([tInd_RFMin tInd_RFMax])) ].' ; 
       
    % Gaussian: Orientation 
    Gabor.Gaussian.Theta(nLR) = Gabor.Sinusoid.Theta(nLR);
    
    % Gaussian: Spread
    Gabor.Gaussian.Spread(:,nLR) =  sqrt(2*pi) * (1/3)* [PatchSize.Ang(1)/2 PatchSize.Ang(2)/2].' ;
    
    % Sinusoid phase is phi_Gabor = arg( FFT(f_opt) )
    Gabor.Sinusoid.Phase(nLR) = RFf.Phase(tInd_FMax) ;
    
end

%% Optimisation

% Constraints
LB0 = [  0 ... 
        -PatchSize.Ang / 2 ... 
        0 ... 
        1./PPA ... 
        1/ max( PatchSize.Ang ) ... 
        0 ...
        0 ... 
     ].';
UB0 = [  Inf ... 
        PatchSize.Ang / 2 ... 
        180 ...
        2*PatchSize.Ang ... 
        (1/2) * max( PPA ) ... 
        180 ... 
        360 ... 
      ].';

Aineq = [];
bineq = [];
          
% Optimise the left and the right images separately. 
for nRun = 1:BestOfNRuns
    
    for nLR = 1:NLR
        
        % Starting point
        P0 = [  Gabor.ScalingFactor(nLR) ... 
                Gabor.Centre(:,nLR).' ...
                Gabor.Gaussian.Theta(:,nLR) ... 
                Gabor.Gaussian.Spread(:,nLR).' ... 
                Gabor.Sinusoid.Frequency(:,nLR) ... 
                Gabor.Sinusoid.Theta(:,nLR) ... 
                Gabor.Sinusoid.Phase(:,nLR) ... 
             ].';
        
        % Redefine the LB and UB based on the general template
        LB = LB0; UB = UB0; 
        
        % Change the frequency bounds. It makes the fitting very specific
        % to Gabors
        LB(7) = P0(7) * 0.25 ; UB(7) = P0(7) * 2 ;
        
        
        % Initialise global optimisation problem
        switch ForceOrientationDependance
            case 0
                Aeq = [] ;
                beq = [] ;                
                
            case 1
                Aeq = [0 0 0 1 0 0 0 -1 0] ;
                beq = 0 ;
        end
        
        % Scaling the error for faster and more accurate calcualtion
        ErrorScalingFactor = 1e-6 ;

        GProb = createOptimProblem( 'fmincon','objective',...
            @(x)fObjFitGabors( x, RFs(:,:,nLR), Mesh, ErrorScalingFactor), ...
            'x0',P0,...
            'Aineq',Aineq,'bineq',bineq,...
            'Aeq',Aeq,'beq',beq ,...
            'lb',LB,'ub',UB');
        
        % Initialise the solver
        GSolver = GlobalSearch('Display','off'); %('UserParallel',true);
        
        % Calculate optima
        [ P_opt_NRuns(:,nLR, nRun) , SS_res_NRuns_By_k(nLR, nRun) ] = run(GSolver,GProb);
        
        SS_res_NRuns(nLR, nRun) = (ErrorScalingFactor * SS_res_NRuns_By_k(nLR, nRun)) ;
        R2_NRuns(nLR, nRun) = 1 -  SS_res_NRuns(nLR, nRun) / SS_tot(nLR) ;
        
    end
    
end

%% Picking the best fits based on R2 amongst the BestOfNRuns runs

% Find the best run
[R2 , NBestRun] = max(R2_NRuns,[],2) ; 

for nLR = 1:NLR
    
    SS_res(nLR) = SS_res_NRuns(nLR, NBestRun(nLR) ) ;

    % Optimal parameters for L/R separately
    P_opt(:,nLR) = squeeze( P_opt_NRuns(:,nLR, NBestRun(nLR) ) ) ;
    
end


%% The optimised parameters, fitting error and R2 goodness-of-fit

FittedGabor.FittingErrors = SS_res ;
FittedGabor.R2 = R2 ;
FittedGabor.R2_NRuns = R2_NRuns ; 
FittedGabor.Parameters_vec = P_opt ;
FittedGabor.Parameters.ScalingFactor        = P_opt(1  ,:) ;
FittedGabor.Parameters.Centre               = P_opt(2:3,:) ;
FittedGabor.Parameters.Gaussian.Theta       = P_opt(4  ,:) ;
FittedGabor.Parameters.Gaussian.Spread      = P_opt(5:6,:) ;
FittedGabor.Parameters.Sinusoid.Frequency   = P_opt(7  ,:) ;
FittedGabor.Parameters.Sinusoid.Theta       = P_opt(8  ,:) ;
FittedGabor.Parameters.Sinusoid.Phase_rel   = P_opt(9  ,:) ;
FittedGabor.Parameters.Sinusoid.Phase_abs   = ...
    (FittedGabor.Parameters.Sinusoid.Phase_rel / 360) .* ...
    (1 ./ FittedGabor.Parameters.Sinusoid.Frequency ) ;

%% Plot the fitted Gabor if required

if ( PlotFittedGabor01 == 1)
    
    % Mesh
    x = [Mesh.XX(:) Mesh.YY(:)].';

    % Calculate the fitted Gabor RFs
    for nLR = 1:NLR
        
        % Gabor
        k = P_opt(1,nLR); 
        c =  P_opt(2:3,nLR); 
        x_centred = x - c*ones(1, size(x,2) ) ; 
        
        G.U = ...
            [cosd(P_opt(4,nLR)) -sind(P_opt(4,nLR)) ; ...
            sind(P_opt(4,nLR))  cosd(P_opt(4,nLR))  ] ; 
        G.S = diag( P_opt(5:6,nLR) ); 
        
        S.f = P_opt(7,nLR); 
        S.V = [ cosd(P_opt(8,nLR)) sind(P_opt(8,nLR)) ].'; 
        S.p = P_opt(9,nLR); 
        
        % RF
        gRF(:,:,nLR) = real( ...
            k * ...
            ( exp( -pi * sum ( (  deg2rad(G.S) \ (G.U.') * deg2rad( x_centred ) ).^2 , 1) ).' ) .* ...
            ( exp( 1i*(2*pi* ( 1/deg2rad(1/S.f) ) * ( S.V.') * deg2rad(x_centred) + deg2rad(S.p) ) ).' )...
                            );
    
    end
    
    % Plot the fitted Gabors
    figure('position',[100 100 1000 600],'color','w');
    
    % Limits for the axes
    MaxSensitivity = max( abs( [RFs(:) ; gRF(:)] ) );
    
    for nLR = 1:NLR
        
        % Actual RFs
        subplot(NLR,2, 1 + (nLR-1)*2 );
        surfc(Mesh.XX,Mesh.YY,RFs(:,:,nLR),...
            'EdgeColor','none');
        view([0 90]);
        colorbar; axis('tight');
        pbaspect([PatchSize.Ang(1) PatchSize.Ang(2) max(PatchSize.Ang)]);
        caxis([-1 1]*MaxSensitivity); colormap(jet);
        zlim([-1 1]*MaxSensitivity) ;
        if (nLR == 1)
            title('Actual RF');
        end
        
        % Fitted Gabor
        subplot(NLR,2, 2 + (nLR-1)*2 );
        surfc(Mesh.XX,Mesh.YY,reshape( gRF(:,:,nLR),size(Mesh.XX)),...
            'EdgeColor','none');
        view([0 90]);
        colorbar; axis('tight');
        pbaspect([PatchSize.Ang(1) PatchSize.Ang(2) max(PatchSize.Ang)]);
        caxis([-1 1]*MaxSensitivity); colormap(jet);
        zlim([-1 1]*MaxSensitivity) ;
        if (nLR == 1)
            title('Fitted Gabor RF');
        end
    end
    
    
end
