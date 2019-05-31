function [X_Center , Y_Center] =fGenerateParameterSubsetPatches(methodROI,NPatches, PatchesSize,ImTrimmedSize,Param )
% The purpose of this function is to generate Random position for NPacthes,
% the output of this function can be used as an argument for fPickPatches.
% % % % % % % % % % % % % % % % % % % 
%       INPUTs
% % % % % % % % % % % % % % % % % % 
% methodROI                 <string> {'Upper', 'Lower', 'Left', 'Right','Eccentricity', 'Manual'} ROI selection mode. All options except 'Manual' are EXPERIMENTAL.
% NPatches                  [1] Number of Generated Patches
% PatchesSize               [2] Size of extracted patches in pixels [X Y]
% ImTrimmedSize.Ang         [2] Size of the Trimmed Image in visual angle
% ImTrimmedSize.Pix         [2] Size of the Trimmed Image in pixels
% Param                     <optional>
% Param.Eccentricity.Min    [1] <Used only if ROI = Eccentricity>; Exemple for select eccentricity between 0 and 2 degreee use :  Param.Eccentricity.Min = 0 & Param.Eccentricity.Max =2;
% Param.Eccentricity.Max    [1] <Used only if ROI = Eccentricity>;
% Param.AvailablePosition   <boolean> [ImTrimmedSize.Pix] Same size than the Trimmed Image, put true for available patches's centers positions
% 
% % % % % % % % % % % % % % % % 
%       OUTPUTs
% % % % % % % % % % % % % % % % 
% X_Center                  [Npatches]
% Y_Center                  [Npatches]
% 
%% Convert position in visual angle cartesian and polar coordinates

PPA  = ImTrimmedSize.Pix ./ ImTrimmedSize.Ang; 
HalfPatchAng = PatchesSize./(2*PPA) ;

% Compute eccentricity for each pixel ( note that this excludes pixels near edges ... )
[XX, YY] = meshgrid( linspace( - ImTrimmedSize.Ang(1)/2 + HalfPatchAng(1), ImTrimmedSize.Ang(1)/2 - HalfPatchAng(1), ImTrimmedSize.Pix(1) - PatchesSize(1) ) , linspace( - ImTrimmedSize.Ang(2)/2 + HalfPatchAng(2), ImTrimmedSize.Ang(2)/2 - HalfPatchAng(2), ImTrimmedSize.Pix(2) - PatchesSize(2)  ) ) ;
YY = flipud(YY);

XX = XX(:); 
YY = YY(:);

Eccentricity  = (XX .^2 + YY .^2).^(1/2);

%% Pick centers of Patches from the subset of available positions

switch lower(methodROI)
    case lower('Upper')
      X_Center = round(randi([ ceil(PatchesSize(1)./2)  ceil(ImTrimmedSize.Pix(1)-PatchesSize(1)./2)], NPatches ,1 ));
      Y_Center = round(randi([ceil(ImTrimmedSize.Pix(2)/2) ImTrimmedSize.Pix(2)], NPatches ,1 ));
      
    case lower('Lower')
      X_Center = round(randi([ ceil(PatchesSize(1)./2)  ceil(ImTrimmedSize.Pix(1)-PatchesSize(1)./2)], NPatches ,1 ));          
      Y_Center = round(randi([0 ceil(ImTrimmedSize.Pix(2)/2)], NPatches ,1 ));
      
    case lower('Left')
      X_Center = round(randi([ceil(PatchesSize(1)./2) ceil(ImTrimmedSize.Pix(1)/2)], NPatches ,1 ));
      Y_Center = round(randi([0  ImTrimmedSize.Pix(2)], NPatches ,1 ));
      
    case lower('Right')
      X_Center = round(randi([ceil(ImTrimmedSize.Pix(1)/2)  ceil(ImTrimmedSize.Pix(1)-PatchesSize(1)./2)], NPatches ,1 ));
      Y_Center = round(randi([0  ImTrimmedSize.Pix(2)], NPatches ,1 ));
      
    case lower('Eccentricity')
      if( exist('Param', 'var') == 0)
          disp(' Error from fSelectPatches: You used "Eccentricity" mode without the associated "Param" ');      
      end
      if (~ isfield(Param, 'Eccentricity')) 
          disp(' Error from fSelectPatches: You used "Eccentricity" mode but without "Param.Eccentricity" ');  
      end
      if ( ( ~ isfield(Param.Eccentricity, 'Min') ) || ( ~ isfield(Param.Eccentricity, 'Max') ) ) 
          disp(' Error from fSelectPatches: You used "Eccentricity" mode without "Param.Eccentricity.Min" & "Param.Eccentricity.Max" ');  
      end

      idxAvailable = find(Eccentricity >= Param.Eccentricity.Min & Eccentricity <= Param.Eccentricity.Max) ;     
      idx = idxAvailable( randi(length(idxAvailable), NPatches,1) ) ; % Pick NPatches random positions from available positions

      [X_Center , Y_Center] = ind2sub([(ImTrimmedSize.Pix(1) - PatchesSize(1)) (ImTrimmedSize.Pix(2) - PatchesSize(2))], idx) ;
      X_Center = X_Center + floor( PatchesSize(1)/2 );
      Y_Center = Y_Center + floor( PatchesSize(2)/2 );
      
    case lower('Manual')
      if( exist('Param', 'var') == 0)
          disp(' Error from fSelectPatches: You used "Manual" mode but without the associated "Param" ');
      end
      if (~ isfield(Param, 'AvailablePosition')) 
          disp(' Error from fSelectPatches: You used "Manual" mode without "Param.AvailablePosition" ');  
      end
      if( (size(Param.AvailablePosition, 2) ~= ImTrimmedSize.Pix(1) ) || (size(Param.AvailablePosition, 1) ~= ImTrimmedSize.Pix(2) ))
          disp(' Error from fSelectPatches: size of "Param.AvailablePosition" should be the same as "ImTrimmedSize.Pix" ') ; 
      end
      
     idxAvailable = find(Param.AvailablePosition);
     idx = idxAvailable( randi(length(idxAvailable), NPatches,1) ) ;
     [X_Center , Y_Center] = ind2sub([ImTrimmedSize.Pix(2) ImTrimmedSize.Pix(1)], idx) ;
      
    otherwise   
    disp('Error from fSelectPatches: The ROI selection mode is invalid. See help for available modes.');
    
end
