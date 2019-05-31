function OutputWeights = fUpdateWeights(InputWeights, Masks, LTP, LTD )
% This function updates the weights of the layer.
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % 
%       INPUTs
% 
% InputWeights               [NSynapses x NNeurones] Synapse weight, specific to each neurons and each inputs 
% 
% Masks.LTP                  [NSynapses x NNeurones] {0,1} A boolean mask describing which synapses must be potentiated 
% Masks.LTD                  [NSynapses x NNeurones] {0,1} A boolean mask describing which synapses must be depressed 
% 
% LTP.Rate                   [NNeurones x NNeurones] LTP rates for each neurone
% LTP.Mu                     [1] Non linearity exponent. It is 0 for purely additive and 1 for purely multiplicative potentiation. 
% LTP.Bounds                 <string> {'soft','hard'} Bounding method. A soft bound uses x*(1-x) as a bound for the potentiation whereas a hard bound strictly implements a (0,1) range for the final weights.
% 
% LTD.Rate                   [NNeurones x NNeurones] LTD rates for each neurone
% LTD.Mu                     [1] Non linearity exponent. It is 0 for purely additive and 1 for purely multiplicative depression. 
% LTD.Bounds                 <string> {'soft','hard'} Bounding method. A soft bound uses x*(1-x) as a bound for depression whereas a hard bound strictly implements a (0,1) range for the final weights.
% 
% % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % 
%       OUTPUTs
% 
% OutputWeights                         [NSynapses x NNeurones] Final synapse weights 
% 

%% Extract relevant variables

% Weights
Ws_vec = InputWeights ; 

%% Calculate potentiation changes, with proper bounds

if LTP.Mu ==0
    dLTP =  Masks.LTP  * LTP.Rate ;

else
    dLTP =  (( (Masks.LTP  .*(1 - Ws_vec) ).^LTP.Mu ) *  LTP.Rate ) ; % diag(

end

% Final weights
switch lower( LTP.Bounds )
    case lower('soft')
        Ws_LTP = Ws_vec + dLTP .* Ws_vec .* (1 - Ws_vec) ;
    case lower('hard')
        Ws_LTP = Ws_vec + dLTP ;
        Ws_LTP(Ws_LTP >= 1) = 1 - realmin ; 
        Ws_LTP(Ws_LTP <= 0) = 0 + realmin ;
end

%% Calculate depression changes

if LTD.Mu ==0
    dLTD =  Masks.LTD  * LTD.Rate ;

else
    dLTD =  ((Masks.LTD  .* Ws_vec ) .^LTD.Mu ) * LTD.Rate ;

end

% Final weights
switch lower( LTD.Bounds )
    case lower('soft')
        Ws_LTD = Ws_vec + dLTD .* Ws_vec .* (1 - Ws_vec) ;
    case lower('hard')
        Ws_LTD = Ws_vec + dLTD ;
        Ws_LTD(Ws_LTD >= 1) = 1 - realmin ;
        Ws_LTD(Ws_LTD <= 0) = 0 + realmin ;
end

%% Calculate final weights
Masks.NoFire = ~(Masks.LTP + Masks.LTD);

OutputWeights = Ws_LTP .* Masks.LTP + Ws_LTD .* Masks.LTD + InputWeights .* Masks.NoFire ;
