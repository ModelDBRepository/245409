function Patches = fPickPatches( TrimmedFilteredImage , PatchesSize , X_Center, Y_Center )
% This function extracts patches from TrimmedFilteredImage, with the given
% centers (X_Center Y_Center ) and the given size (PatchesSize). It works
% best with outputs from fGenerateParameterSubsetPatches.
%
% % % % % % % % % % % % % % % % % % % 
%       INPUTs
% % % % % % % % % % % % % % % % % % 
% TrimmedFilteredImage      [height x width x [] ] The filtered image from which patch should be extracted
% PatchesSize               [2] Size of extracted patches in pixels [X Y]
% X_Center                  [NPatches]
% Y_Center                  [NPatches]
% % % % % % % % % % % % % % % % 
%       OUTPUTs
% % % % % % % % % % % % % % % % 
% Patches                   [ PatchesSize x [] x NPatches]
% 
% 
%% Check inputs

if(length(X_Center) ~= length(Y_Center)) 
    disp('Error from fPickPatches : X_Center and y_center doesnt have the same length ... ') ;
end

%% Resize inputs

dim = size(TrimmedFilteredImage);
TrimmedFilteredImage_Reshape =  reshape(TrimmedFilteredImage, dim(1), dim(2),[] );

%% Extract patches

NPatches = length(X_Center);
Patches = zeros([fliplr(PatchesSize) size(TrimmedFilteredImage_Reshape,3) NPatches]);
for i = 1:NPatches
    
    X_ = X_Center(i); 
    Y_ = Y_Center(i); 

    for d = 1:size(TrimmedFilteredImage_Reshape,3)
        Patches(:,:, d ,i) = TrimmedFilteredImage_Reshape( (X_ - floor((PatchesSize(2) - 1) ./ 2) ):(X_ + ceil((PatchesSize(2) - 1) ./ 2) ) , (Y_- floor((PatchesSize(1) - 1) ./ 2) ):(Y_ + ceil((PatchesSize(1) - 1) ./ 2) ), d );
    end

end


%% Resize outputs

Patches = reshape(Patches, [PatchesSize  dim(3:end) NPatches] );
