clc
clear all 
close all

Task='task';
Group='SN';
NumbROIs=9;
ModelNumb='final';

%add relevant paths
addpath(fullfile(pwd,'/utils'));
% Load task condition 
load('./data/stimulation_condition.mat');
% Load timeseries data 
load('./data/timeseries_Rat_AI_Net_SN_9ROIs.mat');
% Load BSDS model
load('./model/BSDS_model_Rat_AI_Net_SN_9ROIs.mat');
% Load covariance matrix
load('./model/BSDS_covMtx_Rat_AI_Net_SN_9ROIs.mat');

%% Data parameters
TR=0.1; %10Hz
num.ROI=NumbROIs;
num.Subj=length(model.temporal_evolution_of_states); % number of subjects.
num.Vol=length(model.temporal_evolution_of_states{1,1}); % length of timeseries
num.State=length(unique(cell2mat(model.temporal_evolution_of_states)));
num.Run=1;
ROI_names={'PrL' 'AI', 'CG', 'RSC-2.9mm', 'RSC-3.9mm', 'RSC-4.9mm', 'RSC-5.9mm', 'RSC-6.9mm', 'RSC-7.8mm'};


%% Check timeseries data of AI
for subj=1:num.Subj
    AI_ts(subj,:)=data{1,subj}(2,:);
end


%% Reshape group stats (i.e., temporal evolution of state, posterior probability)
% reshape temporal evolution of state into [subject x run] x [time] & identify dominant states
grpTempEvol=reshape(cell2mat(model.temporal_evolution_of_states),[num.Vol,num.Subj*num.Run])';
dominantStates = unique(grpTempEvol(:));
num.State=length(unique(cell2mat(model.temporal_evolution_of_states)));

% relabel states for ease of computation
grpTempEvol_relabel = zeros(size(grpTempEvol));
for relabel = 1:num.State
    grpTempEvol_relabel(grpTempEvol == dominantStates(relabel)) = relabel;
end

% reshape subjects' posterior probability(pp)
subjStatePP = zeros(num.Subj*num.Run, num.Vol, num.State); % each subject's posterior prob.
idx = 1;
for subj =1:num.Subj
    for run = 1:num.Run
        subjStatePP(idx,:,:) = model.posterior_probabilities{1,subj}((num.Vol*run-(num.Vol-1):num.Vol*run), dominantStates);
        idx = idx+1;
    end
end
grpStatePP=squeeze(mean(subjStatePP,1));% compute group posterior probability by averaging subjects posterior prob.


%% Occupancy rates in OFF & ON periods
ON_periods=find(stimulation_cond==1);
OFF_periods=find(stimulation_cond==0);

OR_ON=[]; OR_OFF=[]; OR_Trans=[];
MLT_ON=[]; MLT_OFF=[]; MLT_Trans=[];

for subj=1:num.Subj
    tempEvol_ON=grpTempEvol_relabel(subj,ON_periods);
    [OR_ON(subj,:),MLT_ON(subj,:),~]=summary_stats_fast(tempEvol_ON, 1:num.State);
    tempEvol_OFF=grpTempEvol_relabel(subj,OFF_periods);
    [OR_OFF(subj,:),MLT_OFF(subj,:),~]=summary_stats_fast(tempEvol_OFF, 1:num.State);
end


%% Compute State switching matrix
subjStatSwitchProb = zeros(num.Subj, num.State, num.State); % each subject's state switching matrix
for subj=1:num.Subj
    % y-axis:state(t), x-axis:state(t+1)
    SwitchingMatrix = full(sparse(grpTempEvol_relabel(subj, 1:end-1),grpTempEvol_relabel(subj, 2:end),1,num.State,num.State));
    sum_SwitchingMatrix = sum(SwitchingMatrix,2);
    % Normalize the transition matrix
    for state=1:num.State
        SwitchingMatrix(state, :) = SwitchingMatrix(state,:)/sum_SwitchingMatrix(state);
    end
    subjStatSwitchProb(subj,:,:) = SwitchingMatrix;
end
grpStatSwitchProb = squeeze(nanmean(subjStatSwitchProb(1:num.Subj,:,:)));


%% %%%%%%% Plot figures %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Colorstrings = {'#2F00FF','#FF0A00', '#FFDA0B', '#617296', '#B5B5B5'};
% Colormap=[255 10 9; 255 218 11; 47 0 255; 97 114 150; 181 181 181]./255;
Colormap=[47 0 255; 255 10 9; 255 218 11;  97 114 150; 181 181 181]./255;
    
%% Plot temporal response of AI to optogenetic stimulation
figure;
subplot(2,1,1); plot(stimulation_cond);
xlim([0 length(stimulation_cond)]); title('Stimulation protocol');
subplot(2,1,2); plot(mean(AI_ts));
xlim([0 length(stimulation_cond)]); title('Time course of AI averaged across 9 Chronos rats');


%% reorder the states
OR_ON_re=[]; OR_OFF_re=[];
grpStatePP_re=[];
grpStatSwitchProb_re=[];
grpTempEvol_relabel_re=zeros(size(grpTempEvol_relabel));

for i=1:num.State
    if i==1; i_re=2;        elseif i==2; i_re=4;
    elseif i==3; i_re=1;    elseif i==4; i_re=5;
    elseif i==5; i_re=3;
    end    
    grpStatePP_re(:, i_re) = grpStatePP(:,i);
    OR_ON_re(:,i_re) = OR_ON(:,i);
    OR_OFF_re(:,i_re) = OR_OFF(:,i);    
    
    idx=find(grpTempEvol_relabel==i);
    grpTempEvol_relabel_re(idx) = i_re;

    for j=1:num.State
        if j==1; j_re=2;        elseif j==2; j_re=4;
        elseif j==3; j_re=1;    elseif j==4; j_re=5;
        elseif j==5; j_re=3;
        end        
    	grpStatSwitchProb_re(i_re,j_re)=grpStatSwitchProb(i,j);
    end
end


%% Plot task design & time-varying posterior probability & temporal evolution
figure;
subplot(3,1,1); plot(stimulation_cond);
xlim([0, num.Vol]); title('Task design');

subplot(3,1,2);
plot(grpStatePP);
for i = 1:num.State
    hold on; plot(grpStatePP_re(:,i), 'Color', Colorstrings{i}, 'LineWidth', 2);
end
ylim([0,1]); xlim([0, num.Vol]); xlabel('Time(s)', 'FontSize', 12); ylabel('posterior probability', 'FontSize', 12);
title('Posterior probability of latent brain states');

subplot(3,1,3);
hAxes = gca; imagesc(hAxes, grpTempEvol_relabel_re);
colormap( hAxes, Colormap); title('Temporal Evolution of latent brain states');



%% Plot occupancy rates during ON and OFF stimulation blocks
figure;
subplot(1,2,1);
b=bar(1:num.State, mean(OR_OFF_re)*100, 'FaceColor', 'flat');
for k=1:num.State
    b.CData(k,:) = Colormap(k,:);
end
hold on; errorbar(1:num.State, mean(OR_OFF_re).*100, (std(OR_OFF_re)/sqrt(length(OR_OFF_re))).*100, 'LineStyle', 'none', 'Color', 'k');  % plot standard error.
set(gca, 'xticklabel', {'S1', 'S2', 'S3', 'S4', 'S5'});
ylim([0 100]);
ylabel('Percentage');
title('Occupancy rate during OFF condition');

subplot(1,2,2);
b=bar(1:num.State, mean(OR_ON_re)*100, 'FaceColor', 'flat');
for k=1:num.State
    b.CData(k,:) = Colormap(k,:);
end
hold on; errorbar(1:num.State, mean(OR_ON_re).*100, (std(OR_ON_re)/sqrt(length(OR_ON_re))).*100, 'LineStyle', 'none', 'Color', 'k'); % plot standard error.
set(gca, 'xticklabel', {'S1', 'S2', 'S3', 'S4', 'S5'});
ylim([0 100]);
ylabel('Percentage');
title('Occupancy rate during ON condition');


%% Plot state switching matrix
figure;
imagesc(grpStatSwitchProb);
ylabel('state(t)'); xlabel('state(t+1)');
title('State transition probability');
caxis([0 1]);
textStrings = num2str(grpStatSwitchProb_re(:), '%0.2f'); % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
[x, y] = meshgrid(1:num.State);  % Create x and y coordinates for the strings
hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the strings
                'HorizontalAlignment', 'center');
midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
% Choose white or black for the text color of the strings so they can be easily seen over the background color
textColors = repmat(grpStatSwitchProb_re(:) > midValue, 1, 3);  
set(hStrings, {'Color'}, num2cell(textColors, 2));  % Change the text colors
set(gca, 'XTick', 1:num.State, ...                             % Change the axes tick marks
         'XTickLabel', {'S1', 'S2', 'S3', 'S4', 'S5'}, ...  %   and tick labels
         'YTick', 1:num.State, ...
         'YTickLabel', {'S1', 'S2', 'S3', 'S4', 'S5'}, ...
         'TickLength', [0 0]);
mycolormap = customcolormap([0 0.5 1], {'#9d0142','#f66e45','#ffffbb'});
colorbar('southoutside');
colormap(mycolormap);



