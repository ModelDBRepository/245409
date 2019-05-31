function [Latencies, Diagnostic] = fIntencity2Latency(Intensity, Method, Param)
% This function converts Intensity to spike latency.
% 
% Inputs:
% Intensity     <matrix> The intensity outputs of somes filters. They are
%                        converted to latencies with no change in the size 
%                        of the matrix/array.
% Method        <string> Method of conversion (Only '1/x' is currently 
%                        implemented) and tested. 
% Params        <struct> This structure contains parameters for the method.
% 
% Outputs:
% 
% Latencies                 <matrix> An array/matrix of latencies the same
%                                    size as the intensity array.  
% Diagnostic.Error          <string> A string describing the error.
% Diagnostic.TimeProcess    <scalar> Processing time.

%% Start Clock
Tstart =  tic;

%% Default Method
if nargin < 2
   Method = '1/x'; 
end

%% Implement the intensity to latency conversion

switch lower(Method)
    
    case lower('1/x')
        % Infer 'latencies', with higher output denoting a lower latency.
        % Note that intensity equal to 0 give infinite latency
        Latencies =  1./Intensity;
        
    case lower('ExponantialDecay')
        if nargin < 3
            disp('Error from fIntencity2Latency: you need to supply "Param.Tau" when you use the "ExponantialDecay" method ... ');
        end
        Latencies =  exp( - Intensity /Param.Tau );
        Latencies(Latencies ==1) = Inf;
                
    otherwise
        disp(['ERROR! No mode called ' Method ' exists.']);
        Diagnostic.Error = [ ' "' Method '" is an invalid value for argument "Method" '];
        
end

Diagnostic.TimeProcess = toc(Tstart);