function fplotUpdateWeight( Ws_vecs_TimeCourse, Binocular, Labels )
% This function plots the weights of the synapses for each neuron of a given layer. The purpose of the plot is to see the update and the convergence of the weight.
%
% INPUTs:
% Ws_vecs_TimeCourse      [Synapses x Neurons x Times] The time course of the Weights
% Binocular               <boolean>  The layer is binocular ?
% Labels                  [Neurones x 1] A numerical label for each of the neurones to be plotted
%
% OUTPUTs
% [void]

%% Extract parameters
W0_vec = squeeze(Ws_vecs_TimeCourse(:,:,1));
NFeatures = size(Ws_vecs_TimeCourse,1);
TrainingSet.Nframe = size(Ws_vecs_TimeCourse,3);
NNeurones = size(Ws_vecs_TimeCourse,2);


%%
if Binocular == 1
    
    % Number of snapshots of how the fields developed.
    NSnapshots = TrainingSet.Nframe + 1 ; 
    FilmDuration = 30; %seconds
    
    fps = NSnapshots/FilmDuration ;
    
    % Indices of the frame
    tIndz = ceil( linspace(1, TrainingSet.Nframe, NSnapshots-1) );
    
    for nnNeurone = 1:5:NNeurones
        % Plot the initial weights, five at a time
        hFiggins( (4+nnNeurone)/5 ) = figure('Position',[10 100 1400 700],'Visible','On');
        
        for nnnNeurone = 1:5
            
            nNeurone = nnnNeurone + nnNeurone - 1;
            if (nNeurone <= NNeurones)
                
                subplot(1,5,nnnNeurone );
                
                plot(1:NFeatures/2,...
                    W0_vec((1:NFeatures/2),nNeurone),'bo',...
                    1:NFeatures/2,...
                    W0_vec((1:NFeatures/2) + (NFeatures/2),nNeurone),'rx');
                hold on;
                
                tkde_L = ksdensity(W0_vec((1:NFeatures/2),nNeurone),linspace(0,1,10));
                tkde_L = (NFeatures/2)* tkde_L / sum(tkde_L(:));
                tkde_R = ksdensity(W0_vec((1:NFeatures/2) + (NFeatures/2),nNeurone),linspace(0,1,10));
                tkde_R = (NFeatures/2)*tkde_R / sum(tkde_R(:));
                plot(tkde_L,linspace(0,1,10),'b--',...
                    tkde_R,linspace(0,1,10),'r--',...
                    'LineWidth',2);
                hold on;
                
                L1s = sum(W0_vec((1:NFeatures/2),nNeurone) >= 0.5);
                R1s = sum(W0_vec((1:NFeatures/2) + (NFeatures/2),nNeurone) >= 0.5);
                legend({[num2str(L1s) ' > 0.5, ' num2str(NFeatures/2 - L1s) ' < 0.5' ],...
                    [num2str(R1s) ' > 0.5, ' num2str(NFeatures/2 - R1s) ' < 0.5' ]},...
                    'Location','SouthOutside');
                
                title(['Neurone: ' num2str(Labels(nNeurone)) '; t: 0']);
                axis('tight'); ylim([0 1]);
                hold off;
                
            end
        end
    end
    
    pause(0.001/fps);
    % tFilm(1) = getframe(gcf);
    
    for ntt = 1:numel(tIndz)
        tt = tIndz(ntt);
        
        for nnNeurone = 1:5:NNeurones
            set(0,'CurrentFigure', hFiggins( (4+nnNeurone)/5 ) ) ;
            for nnnNeurone = 1:5
                nNeurone = nnnNeurone + nnNeurone - 1;
                
                if (nNeurone <= NNeurones)
                    
                    subplot(1,5,nnnNeurone );
                    
                    tWs_L_t = Ws_vecs_TimeCourse((1:NFeatures/2),nNeurone,tt) ;
                    tWs_R_t = Ws_vecs_TimeCourse((1:NFeatures/2) + (NFeatures/2),nNeurone,tt) ;
                    plot(1:NFeatures/2,tWs_L_t,'bo',...
                        1:NFeatures/2,tWs_R_t,'rx');
                    hold on;
                    
                    tkde_L = ksdensity(tWs_L_t,linspace(0,1,10));
                    tkde_L = (NFeatures/2)* tkde_L / sum(tkde_L(:));
                    tkde_R = ksdensity(tWs_R_t,linspace(0,1,10));
                    tkde_R = (NFeatures/2)*tkde_R / sum(tkde_R(:));
                    plot(tkde_L,linspace(0,1,10),'b--',...
                        tkde_R,linspace(0,1,10),'r--',...
                        'LineWidth',2);
                    hold on;
                    
                    
                    L1s = sum(Ws_vecs_TimeCourse((1:NFeatures/2),nNeurone,tt) >= 0.5);
                    R1s = sum(Ws_vecs_TimeCourse((1:NFeatures/2) + (NFeatures/2),nNeurone,tt) >= 0.5);
                    legend({[num2str(L1s) ' > 0.5, ' num2str(NFeatures/2 - L1s) ' < 0.5' ],...
                        [num2str(R1s) ' > 0.5, ' num2str(NFeatures/2 - R1s) ' < 0.5' ]},...
                        'Location','SouthOutside');
                    
                    title(['Neurone: ' num2str(Labels(nNeurone)) '; t: ' num2str(tt)]);
                    axis('tight'); ylim([0 1]);
                    hold off;
                    
                end
            end
        end
        % drawnow;
        pause(0.001/fps);
    end
    
else
    % Number of snapshots of how the fields developed.
    NSnapshots = TrainingSet.Nframe/500+1;
    FilmDuration = 15; %seconds
    
    fps = NSnapshots/FilmDuration ;
    
    % Indices of the frame
    tIndz = ceil( linspace(1,TrainingSet.Nframe,NSnapshots-1) );
    
    for nnNeurone = 1:5:NNeurones
        hFiggins( (4+nnNeurone)/5 ) = figure('Position',[10 100 1400 700],'Visible','On');
        
        for nnnNeurone = 1:5%NNeurones%
            
            nNeurone = nnnNeurone + nnNeurone - 1;
            if (nNeurone <= NNeurones)
                
                subplot(1,5,nnnNeurone );
                
                plot(1:NFeatures, W0_vec(1:NFeatures,nNeurone),'bo');
                
                hold on;
                
                tkde_L = ksdensity(W0_vec((1:NFeatures),nNeurone),linspace(0,1,10));
                tkde_L = (NFeatures)* tkde_L / sum(tkde_L(:));
                
                plot(tkde_L,linspace(0,1,10),'b--','LineWidth',2);
                hold on;
                
                L1s = sum(W0_vec((1:NFeatures),nNeurone) >= 0.5);
                legend({[num2str(L1s) ' > 0.5, ' num2str(NFeatures - L1s) ' < 0.5' ]},'Location','SouthOutside');
                title(['Neurone: ' num2str(nNeurone) '; t: 0']);
                axis('tight'); ylim([0 1]);
                hold off;
                
            end
        end
    end
    
    pause(0.001/fps);
    
    for ntt = 1:numel(tIndz)
        tt = tIndz(ntt);
        
        for nnNeurone = 1:5:NNeurones
            set(0,'CurrentFigure', hFiggins( (4+nnNeurone)/5 ) ) ;
            for nnnNeurone = 1:5
                nNeurone = nnnNeurone + nnNeurone - 1;
                
                if (nNeurone <= NNeurones)
                    
                    subplot(1,5,nnnNeurone );
                    
                    tWs_L_t = Ws_vecs_TimeCourse((1:NFeatures),nNeurone,tt) ;
                    plot(1:NFeatures,tWs_L_t,'bo');
                    hold on;
                    
                    tkde_L = ksdensity(tWs_L_t,linspace(0,1,10));
                    tkde_L = (NFeatures)* tkde_L / sum(tkde_L(:));
                    plot(tkde_L,linspace(0,1,10),'b--','LineWidth',2);
                    hold on;
                    
                    
                    L1s = sum(Ws_vecs_TimeCourse((1:NFeatures),nNeurone,tt) >= 0.5);
                    legend({[num2str(L1s) ' > 0.5, ' num2str(NFeatures - L1s) ' < 0.5' ]},'Location','SouthOutside');
                    
                    title(['Neurone: ' num2str(nNeurone) '; t: ' num2str(tt)]);
                    axis('tight'); ylim([0 1]);
                    hold off;
                    
                end
            end
        end
        pause(0.001/fps);
    end
    
    
end