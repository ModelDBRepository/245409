function Horizontal_DTC = fHorizontalDTC_Correlation( BinocularRF, PPA, PatchSize )
% Function to calculate the DTCs based on BC.
% 
% Inputs:
% BinocularRF: Binocular RFs
% PPA: Pixels per angle
% PatchSize: Size of the patches
% 
% Outputs: 
% Horizontal_DTC: DTCs.


% Extract the left and right RFs
RF.L = squeeze( BinocularRF(:,:,1) ); 
RF.R = squeeze( BinocularRF(:,:,2) ); 

% Calculate cross correlation (overkill, but much easier/quicker to code...)
XCorr_raw = xcorr2( RF.R, RF.L); 

% DTC
Horizontal_DTC.DTC_raw = XCorr_raw( size(RF.R ,1), :).' ; 

% Disparity range
Horizontal_DTC.DisparityRange = 1/PPA(1) * ...
    ( ( 1:numel(Horizontal_DTC.DTC_raw) ) - PatchSize.Pix(1) ) ;

% Fit a Gabor to it
P0 = [0.9*max(abs(Horizontal_DTC.DTC_raw)) 0 0.5 1 randi([1 360]) 0] .' ;
LB = [0   -PatchSize.Ang(1)/2  1/PPA(1)             0.1     0  -max(abs(Horizontal_DTC.DTC_raw))] ;
UB = [Inf  PatchSize.Ang(1)/2  2*PatchSize.Ang(1)  10     360   max(abs(Horizontal_DTC.DTC_raw))] ;
ErrorScaling = 1e16 ; % 1e12
Horizontal_DTC.GaborFit_vec = fmincon( @(x)fObjFit1DGabors( x, ...
                                                 Horizontal_DTC.DTC_raw, ...
                                                 Horizontal_DTC.DisparityRange,...
                                                 ErrorScaling),...
                              P0 , [],[],[],[],...
                              LB, UB) ;
                          
% Fitted gabor points
x = Horizontal_DTC.DisparityRange ;
k = Horizontal_DTC.GaborFit_vec(1) ;
c =  Horizontal_DTC.GaborFit_vec(2) ; 
x_centred = x - c ; 

G.S = Horizontal_DTC.GaborFit_vec(3) ;

S.f = Horizontal_DTC.GaborFit_vec(4); 
S.p = Horizontal_DTC.GaborFit_vec(5); 

Baseline = Horizontal_DTC.GaborFit_vec(6); 

Horizontal_DTC.GaborFit = Baseline + real( k * ...
      ( exp( -pi * ( deg2rad( x_centred ) / deg2rad( G.S ) ).^2  ).' ) .* ...
      ( exp( 1i*( 2*pi* ( 1/deg2rad(1/S.f) ) * deg2rad(x_centred) + deg2rad(S.p) ) ).' ) );