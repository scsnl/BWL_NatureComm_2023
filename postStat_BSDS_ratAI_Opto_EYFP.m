clc
clear all 
close all

Task='task';
Group='EYFP';
NumbROIs=9;
ModelNumb='final';

%add relevant paths
addpath(fullfile(pwd,'/utils'));

% Load task condition 
load(fullfile(pwd, '/data/', '/stimulation_condition.mat'));

% Load model
load(fullfile(pwd, '/model/', sprintf('BSDS_model_ratAI_Opto_%s.mat', Group)));



%% Data parameters
TR = 1; 
num.ROI = NumbROIs;
num.Subj = length(model.temporal_evolution_of_states); % number of subjects.
num.Vol = length(model.temporal_evolution_of_states{1,1}); % length of timeseries
num.State = length(unique(cell2mat(model.temporal_evolution_of_states)));
num.Run = 1;
num.Trials = 8; %number of ON stimulation blocks
ROI_names = {'PrL' 'AI', 'CG', 'RSC-2.9mm', 'RSC-3.9mm', 'RSC-4.9mm', 'RSC-5.9mm', 'RSC-6.9mm', 'RSC-7.8mm'};
ROI_names_excld_AI = {'PrL', 'CG', 'RSC-2.9mm', 'RSC-3.9mm', 'RSC-4.9mm', 'RSC-5.9mm', 'RSC-6.9mm', 'RSC-7.8mm'};
Colorstrings = {'#454246','#e0a96d', '#ddc3a5'};
Colormap=[0.2706 0.2588 0.2745;0.8784 0.6627 0.4275; 0.8667 0.7647 0.6471];
reorder=[3,1,5,4,2]; %reorder the states for the visualization purpose.
    
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


%% Figure 3B: Plot task design & time-varying posterior probability & temporal evolution

close all;
figure;
subplot(3,1,1);
plot(stimulation_cond);
xlim([0, num.Vol]);
title('Task design');

subplot(3,1,2);
plot(grpStatePP);
for i = 1:num.State
    hold on;
    plot(grpStatePP(:,i), 'Color', Colorstrings{i}, 'LineWidth', 2);
end
ylim([0,1]); xlim([0, num.Vol]); xlabel('Time(s)', 'FontSize', 12); ylabel('posterior probability', 'FontSize', 12);
title('Posterior probability of latent brain states');

subplot(3,1,3);
hAxes = gca;
imagesc(hAxes, grpTempEvol_relabel);
colormap( hAxes, Colormap);
title('Temporal Evolution of latent brain states');

%% Figure 3E: Plot occupancy rates
figure;
subplot(1,2,1);
b=bar(1:num.State, mean(OR_OFF)*100, 'FaceColor', 'flat');
for k=1:num.State
    b.CData(k,:) = Colormap(k,:);
end
hold on; errorbar(1:num.State, mean(OR_OFF).*100, (std(OR_OFF)/sqrt(length(OR_OFF))).*100, 'LineStyle', 'none', 'Color', 'k');  % plot standard error.
set(gca, 'xticklabel', {'S1', 'S2', 'S3', 'S4'});
ylim([0 100]);
ylabel('Percentage');
title('Occupancy rate during OFF condition');

subplot(1,2,2);
b=bar(1:num.State, mean(OR_ON)*100, 'FaceColor', 'flat');
for k=1:num.State
    b.CData(k,:) = Colormap(k,:);
end
hold on; errorbar(1:num.State, mean(OR_ON).*100, (std(OR_ON)/sqrt(length(OR_ON))).*100, 'LineStyle', 'none', 'Color', 'k'); % plot standard error.
set(gca, 'xticklabel', {'S1', 'S2', 'S3', 'S4'});
ylim([0 100]);
ylabel('Percentage');
title('Occupancy rate during ON condition');


