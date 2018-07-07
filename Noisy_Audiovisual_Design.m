clc
clear all
close all
path = pwd;

% stimulus names for filename output
stim_con = {'A','V','AV','AN','AVN'};

% read optseq design files
design_files = dir('*.par');
design_files = {design_files(:).name}';

numruns = 4; % number of runs
numstim = 40; % number of stimuli per run

for i = 1:numel(design_files)
    
    [stimonset,stimtype,stimlength,~,stimname] =  textread(design_files{i},'%f%n%f%f%s','delimiter',' ');
    designFile = [num2cell(stimonset),num2cell(stimtype),num2cell(stimlength),stimname];
    
    for j = 1:numel(stim_con)
        
        onset.(stim_con{j}) = zeros(numel(find(stimtype==j)),1);
        duration.(stim_con{j}) = zeros(numel(find(stimtype==j)),1);
        count.(stim_con{j}) = 1;
        
    end
    
    null_onset = zeros(numel(find(stimtype==0)),1);
    null_duration = zeros(numel(find(stimtype==0)),1);
    null_count = 1;
    
    % Create onset vectors
    for j = 1:numel(stimonset)

        if stimtype(j) == 0
            null_onset(null_count) = stimonset(j);
            null_duration(null_count) = stimlength(j);
            null_count = null_count + 1;
        else
            onset.(stim_con{stimtype(j)})(count.(stim_con{stimtype(j)})) = stimonset(j);
            duration.(stim_con{stimtype(j)})(count.(stim_con{stimtype(j)})) = stimlength(j);
            count.(stim_con{stimtype(j)}) = count.(stim_con{stimtype(j)}) + 1;
        end
        
    end
    
    onset.NULL = null_onset;
    duration.NULL = null_duration;
    
    % Save onsets of all designs for conditions separately
    designs_exp.(['design_' num2str(i)]).onset = onset;
    designs_exp.(['design_' num2str(i)]).duration = duration;
    
    designs_exp.(['design_' num2str(i)]).onsets_tot = stimonset';
    designs_exp.(['design_' num2str(i)]).stimduration_tot = stimlength';
    
    designs_exp.(['design_' num2str(i)]).design_tot = [stimtype stimonset]';    
end

clearvars -except designs_exp path stim_con

%% Assign stimuli to conditions
video_files = dir('/Volumes/data/BCM/Experiments/Noisy_Audiovisual_fMRI/stimuli/*.mp4');
video_files = {video_files(:).name}';

% Stimulus info
noise_level = [1,0,1,2,2];
stim_con_noise = find(noise_level==2);
stim_con_clear = find(noise_level==1);
stim_con_silent = find(noise_level==0);

video_files_c1 = video_files(1:32); % Only clear A, AV speech
video_files_n1 = video_files(33:64); % Only noisy A, AV speech

video_files_c2 = video_files(65:72); % Only clear A, AV speech for first repetition
video_files_n2 = video_files(73:end); % Only noisy A, AV speech for first repetition

videos_run_audio = (1:(numel(video_files)/2)-8)';
s = RandStream('mt19937ar','Seed',1);
run_1 = sort(randperm(s,numel(videos_run_audio),numel(videos_run_audio)/2))'; % Stimuli to select for run 1 (from both, clear and noisy, samples)
run_2 = setdiff(videos_run_audio,run_1); % Stimuli to select for run 2 (from both, clear and noisy, samples)

% Select movies to repeat in second run
%remove_rep = sort(randperm(s,numel(run_1),numel(video_files_c2)/2))';
remove_rep = (numel(video_files_c1)/2-3:numel(video_files_c1)/2)';

run_1_mod = run_1;
run_1_mod(remove_rep) = [];

run_2_mod = run_2;
run_2_mod(remove_rep) = [];

% Select movies to use in run 3 for first time
run_3_add = sort(randperm(s,numel(video_files_c2),numel(video_files_c2)/2))';

% Select movies to use in run 4 for first time
run_4_add = setdiff(1:numel(video_files_c2),run_3_add);

clear s

% Assign stimuli first two runs
for i = 1:4
    
    design_c = designs_exp.(['design_' num2str(i)]).design_tot(1,:)';
    n_cond = unique(design_c(1,:));
    n_cond(n_cond==0) = [];
    
    % Select 
    if i == 1
        vid_c = video_files_c1(run_1);
        vid_n = video_files_n1(run_1);
        
        vid_c_a = vid_c(1:2:end);
        vid_c_av = vid_c(2:2:end);
        vid_n_a = vid_n(1:2:end);
        vid_n_av = vid_n(2:2:end);
    elseif i == 2
        vid_c = video_files_c1(run_2);
        vid_n = video_files_n1(run_2);
        
        vid_c_a = vid_c(1:2:end);
        vid_c_av = vid_c(2:2:end);
        vid_n_a = vid_n(1:2:end);
        vid_n_av = vid_n(2:2:end);
    elseif i == 3
        vid_c = [video_files_c1(run_1_mod); video_files_c2(run_3_add)];
        vid_n = [video_files_n1(run_1_mod); video_files_n2(run_3_add)];
        
        vid_c_a = vid_c(2:2:end);
        vid_c_av = vid_c(1:2:end);
        vid_n_a = vid_n(2:2:end);
        vid_n_av = vid_n(1:2:end);
    elseif i == 4
        vid_c = [video_files_c1(run_2_mod); video_files_c2(run_4_add)];
        vid_n = [video_files_n1(run_2_mod); video_files_n2(run_4_add)];
        
        vid_c_a = vid_c(2:2:end);
        vid_c_av = vid_c(1:2:end);
        vid_n_a = vid_n(2:2:end);
        vid_n_av = vid_n(1:2:end);
    end
    
    vid_idx = randi([1 numel(video_files)],1,8);
    vid_v = video_files(vid_idx);
        
    % Assign videos to design
    stimuli_c = cell(numel(design_c),1);
    
    stimuli_c(design_c==1) = vid_c_a; % Assign clear auditory only
    stimuli_c(design_c==2) = vid_v; % Assign visual only
    stimuli_c(design_c==3) = vid_c_av; % Assign clear audiovisual
    stimuli_c(design_c==4) = vid_n_a; % Assign noisy auditory only
    stimuli_c(design_c==5) = vid_n_av; % Assign noisy audiovisual
    stimuli_c(design_c==0) = {'NULL'}; % Assign NULL events
    
    % Assign data type to videos (wav/noisy condition)
    for j = 1:numel(design_c)
        if design_c(j) == 1 % in case of clear auditory only
            stim_c = stimuli_c{j};
            stim_c = strrep(stim_c,'.mp4','.wav');
            stimuli_c{j} = stim_c;
        elseif design_c(j) == 4 % in case of noisy auditory only
            stim_c = stimuli_c{j};
            stim_c = strrep(stim_c,'.mp4','_-8dB.wav');
            stimuli_c{j} = stim_c;
        elseif design_c(j) == 5 % in case of noisy audiovisual
            stim_c = stimuli_c{j};
            stim_c = strrep(stim_c,'.mp4','_-8dB.mp4');
            stimuli_c{j} = stim_c;
        end
    end

    % Current design without NULL as video code for conditions
    design_c_n = design_c;
    design_c_n(design_c_n==0) = [];
    
    designs_exp.(['design_' num2str(i)]).stimuli = stimuli_c;
    designs_exp.(['design_' num2str(i)]).video_code_con = design_c_n;
end

save('Noisy_Audiovisual_fMRI_Design.mat','designs_exp');
