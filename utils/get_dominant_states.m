function [dominant_states, num_state]=get_dominant_states_v2(model, runName)


if strcmp(runName, 'AllRuns')
    num_Runs = 4; % number of runs. Total of 4 (runa, runb, runc, rund)
else
    num_Runs =1;
end
num_Subj = length(model.model.posterior_probabilities); % number of subjects. 
num_Vol = size(model.model.posterior_probabilities{1}, 1)/num_Runs; %length of time series data

grp_state_mtx_by_subj = zeros(num_Subj, num_Vol*num_Runs); %[Number of Subject] x [Volume of each run x Four runs]

for subj=1:num_Subj
    stat_trans_by_subj=[];
	stat_trans_by_subj = model.model.temporal_evolution_of_states{subj};
    grp_state_mtx_by_subj(subj,:)=stat_trans_by_subj;
end

dominant_states = unique(grp_state_mtx_by_subj(:));
num_state = length(dominant_states);