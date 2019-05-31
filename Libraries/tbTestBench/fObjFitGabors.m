function FittingErrorBy1e4 = fObjFitGabors( GaborParameters, RF, Mesh, ErrorScalingFactor)
% This objective function calculates the squared fitting error for a gabor
% given the real receptive field.

%% Extract parameters

% Mesh
x = [Mesh.XX(:) Mesh.YY(:)].';

% Gabor
k = GaborParameters(1); 
c =  GaborParameters(2:3); 
x_centred = x - c*ones(1, size(x,2) ) ; 
G.U = ...
    [cosd(GaborParameters(4)) -sind(GaborParameters(4)) ; ...
     sind(GaborParameters(4))  cosd(GaborParameters(4))  ] ; 
G.S = diag( GaborParameters(5:6) );

S.f = GaborParameters(7); 
S.V = [ cosd(GaborParameters(8)) sind(GaborParameters(8)) ].'; 
S.p = GaborParameters(9); 

%% RF
gRF = k * ...
      ( exp( -pi * sum ( (  deg2rad(G.S) \ (G.U.') * deg2rad( x_centred ) ).^2 , 1) ).' ) .* ...
      ( exp( 1i*( 2*pi* ( 1/deg2rad(1/S.f) ) * ( S.V.') * deg2rad(x_centred) + deg2rad(S.p) ) ).' ) ;
 
%% Calculate the squared fitting error (biased estimator)

tk = ErrorScalingFactor ;
FittingErrorBy1e4 =  mean( ( real(gRF) - RF(:) ).^2 ) / tk ;