function DoG = fCreateDoG( PixelSize, PPA, Sigmas, OnOff)

%% Size: MUST be odd number of pixels. This reduces trimming complexity.

DoG.Size.Pix = PixelSize; % Horiz x Vert
DoG.Size.Ang = DoG.Size.Pix ./ PPA ; % Horiz x Vert

%% Define centre

DoG.Centre.Sigma.Pix = Sigmas(1) ; 
DoG.Centre.Sigma.Ang = DoG.Centre.Sigma.Pix ./ PPA ;
DoG.Centre.Kernel = fspecial('gaussian',fliplr(DoG.Size.Pix),DoG.Centre.Sigma.Pix) ;

%% Define surround

DoG.Surround.Sigma.Pix = Sigmas(2) ; 
DoG.Surround.Sigma.Ang = DoG.Surround.Sigma.Pix ./ PPA ;
DoG.Surround.Kernel = fspecial('gaussian',fliplr(DoG.Size.Pix),DoG.Surround.Sigma.Pix) ;

%% Define the kernel

% Depending on whether it is ON-centre OFF-surround or vice versa. 
switch lower(OnOff)
    case lower('on')
        DoG.Kernel = + DoG.Centre.Kernel - DoG.Surround.Kernel ;
    case lower('off')
        DoG.Kernel = - DoG.Centre.Kernel + DoG.Surround.Kernel ;
end

%% Normalisation: Sum-to-zero
DoG.Kernel = DoG.Kernel - mean( DoG.Kernel(:) ) ;

%% Normalisation: Max output = 1. 

% This normalisation depends on the range of the input values. A typical
% case for normalised luminance is the range [0 1].
Input.Min = 0; Input.Max = 1; 
DoG.Kernel = DoG.Kernel / ...   
    ( ( sum( abs( DoG.Kernel(:) ) ) / 2 ) * (Input.Max - Input.Min) );
