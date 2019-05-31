function FittingError = fObjFit1DGabors( GaborParameters, Activations, Disparities, ErrorScaling )
% This objective function calculates the squared fitting error for a gabor
% given the real receptive field.

%% Extract parameters

% Mesh
x = Disparities;

% Gabor
k = GaborParameters(1) ; 
c =  GaborParameters(2) ;
x_centred = x - c ; 

G.S = GaborParameters(3) ; 

S.f = GaborParameters(4); 
S.p = GaborParameters(5); 

Baseline = GaborParameters(6);

%% Calculate the Gabor RF
gRF = Baseline + ( k * ...
      ( exp( -pi * ( deg2rad( x_centred ) / deg2rad( G.S ) ).^2  ).' ) .* ...
      ( exp( 1i*( 2*pi* ( 1/deg2rad(1/S.f) ) * deg2rad(x_centred) + deg2rad(S.p) ) ).' ) );
 
%% Calculate the squared fitting error (biased estimator)
FittingError =  ErrorScaling * mean( ( real(gRF) - Activations(:) ).^2 ) ;
