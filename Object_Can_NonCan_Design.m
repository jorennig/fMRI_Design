tic

clc
clear all
path = pwd;

% stimulus names for filename output
stimtypes = {'o1_can','o1_noncan','o2_can','o2_noncan'};

% read optseq design files
data_files = dir('*.par');
for i=1:length(data_files)    
    design_files{i} = data_files(i).name;       
end
clear data_files

% number of runs
numruns = numel(design_files);
% number of stimuli per run
numstim = 168;
% how long each stimulus lasts
stimtime = 2;
% max number of stim 
totalstims = numruns*numstim;
% Array containing stim types
stimseq = cell(numruns, numstim);
% Array containing clip names to save to disk
clipfile = cell(numruns,numstim);

% change to directory that contains files
currentPath = pwd;

% Stimulus prefix (from optseq2)
designs = cell(numel(design_files),numstim);
pre_rec = 0;

for i = 1:numel(design_files)
    
    [stimonset,stimtype,stimlength,~,stimname] =  textread(design_files{i},'%f%n%f%f%s','delimiter',' ');
    designFile = [num2cell(stimonset),num2cell(stimtype),num2cell(stimlength),stimname];
    
    % Accounts in seconds for dummy scans in the beginning
    stimonset = stimonset + pre_rec;
    
    o1_can_onset = zeros(numel(find(stimtype==1)),1);
    o1_noncan_onset = zeros(numel(find(stimtype==2)),1);
    o2_can_onset = zeros(numel(find(stimtype==3)),1);
    o2_noncan_onset = zeros(numel(find(stimtype==4)),1);

    null_onset = zeros(numel(find(stimtype==0)),1);
    
    o1_can_duration = zeros(numel(find(stimtype==1)),1);
    o1_noncan_duration = zeros(numel(find(stimtype==2)),1);
    o2_can_duration = zeros(numel(find(stimtype==3)),1);
    o2_noncan_duration = zeros(numel(find(stimtype==4)),1);

    null_duration = zeros(numel(find(stimtype==0)),1);
    
    o1_can_count = 1;
    o1_noncan_count = 1;

    o2_can_count = 1;
    o2_noncan_count = 1;
    
    null_count = 1;
    
    % Create onset vectors
    for j = 1:numel(stimonset)

        if stimtype(j) == 1
            o1_can_onset(o1_can_count) = stimonset(j);
            o1_can_duration(o1_can_count) = stimlength(j);
            o1_can_count = o1_can_count + 1;
        end
        if stimtype(j) == 2
            o1_noncan_onset(o1_noncan_count) = stimonset(j);
            o1_noncan_duration(o1_noncan_count) = stimlength(j);
            o1_noncan_count = o1_noncan_count + 1;
        end
        if stimtype(j) == 3
            o2_can_onset(o2_can_count) = stimonset(j);
            o2_can_duration(o2_can_count) = stimlength(j);
            o2_can_count = o2_can_count + 1;
        end
        if stimtype(j) == 4
            o2_noncan_onset(o2_noncan_count) = stimonset(j);
            o2_noncan_duration(o2_noncan_count) = stimlength(j);
            o2_noncan_count = o2_noncan_count + 1;
        end
        if stimtype(j) == 0
            null_onset(null_count) = stimonset(j);
            null_duration(null_count) = stimlength(j);
            null_count = null_count + 1;
        end
                
    end
    
    % Save onsets of all conditions separately
    designs_exp.(['design_' num2str(i)]).onsets_o1_can = o1_can_onset';
    designs_exp.(['design_' num2str(i)]).onsets_o1_noncan = o1_noncan_onset';
    designs_exp.(['design_' num2str(i)]).onsets_o2_can = o2_can_onset';
    designs_exp.(['design_' num2str(i)]).onsets_o2_noncan = o2_noncan_onset';
    designs_exp.(['design_' num2str(i)]).onsets_null = null_onset';
    
    designs_exp.(['design_' num2str(i)]).duration_o1_can = o1_can_duration';
    designs_exp.(['design_' num2str(i)]).duration_o1_noncan = o1_noncan_duration';
    designs_exp.(['design_' num2str(i)]).duration_o2_can = o2_can_duration';
    designs_exp.(['design_' num2str(i)]).duration_o2_noncan = o2_noncan_duration';
    designs_exp.(['design_' num2str(i)]).duration_null = null_duration';

    % Save onsets for McGurk and Control trials
    in_onset = sort([o1_can_onset; o2_can_onset]);
    sc_onset = sort([o1_noncan_onset; o2_noncan_onset]);

    designs_exp.(['design_' num2str(i)]).onsets_in = in_onset';
    designs_exp.(['design_' num2str(i)]).onsets_sc = sc_onset';
    
    designs_exp.(['design_' num2str(i)]).onsets_exp = stimonset';
    designs_exp.(['design_' num2str(i)]).stimduration_exp = stimlength';
    
    designs_exp.(['design_' num2str(i)]).design = [stimtype stimonset]';    
end

clearvars -except designs_exp path

%% Assign stimuli to conditions
cd('/Users/beauchamplab/Dropbox/_Papers/_New_Projects/Paradigms/fMRI_Object_Can_NonCan/stimuli');

% Read in stimuli
data_files = dir('*C_1.jpg');
for i=1:length(data_files)    
    img_files{i} = data_files(i).name;       
end
img_f.img_files_can_1 = img_files';
clear data_files i

data_files = dir('*N_1.jpg');
for i=1:length(data_files)    
    img_files{i} = data_files(i).name;       
end
img_f.img_files_noncan_1 = img_files';
clear data_files i

data_files = dir('*C_2.jpg');
for i=1:length(data_files)    
    img_files{i} = data_files(i).name;       
end
img_f.img_files_can_2 = img_files';
clear data_files i

data_files = dir('*N_2.jpg');
for i=1:length(data_files)    
    img_files{i} = data_files(i).name;       
end
img_f.img_files_noncan_2 = img_files';
clear data_files i img_files

cd(path);

% Assign stimuli
for i = 1:numel(fieldnames(designs_exp))
    
    fprintf('Create design %d \n',i)
    
    check_r = 0;
    count = 1;
    while check_r == 0

        count = count + 1;
        
        design_c = designs_exp.(['design_' num2str(i)]).design;
        cond = unique(design_c(1,:));
        cond(1) = [];

        for j = 1:numel(cond)

            if j == 1

                a = img_f.img_files_can_1;
                cond_pos = find(design_c(1,:)==cond(j));
                num_el = numel(cond_pos);
                rep = numel(cond_pos)/numel(a);

                code_c = 1:numel(a);
                code_r = repmat(code_c',rep,1);

                check = 0;
                while check == 0

                    code_s = Shuffle(code_r);
                    f = [];
                    for k = 1:length(code_c)
                        pat = findpattern(code_s,ones(2,1)*code_c(k));

                        if isempty(pat)
                            f(k) = 0;
                        else
                            f(k) = 1;
                        end
                    end

                    if sum(f) == 0
                        check = 1;
                    end

                end

                stim_c = a(code_s);
                stimuli_design.(['design_' num2str(i)]).(['c_' num2str(cond(j))]) = stim_c;

            end

            if j == 2

                a = img_f.img_files_noncan_1;
                cond_pos = find(design_c(1,:)==cond(j));
                num_el = numel(cond_pos);
                rep = numel(cond_pos)/numel(a);

                code_c = 1:numel(a);
                code_r = repmat(code_c',rep,1);

                check = 0;
                while check == 0

                    code_s = Shuffle(code_r);
                    f = [];
                    for k = 1:length(code_c)
                        pat = findpattern(code_s,ones(2,1)*code_c(k));

                        if isempty(pat)
                            f(k) = 0;
                        else
                            f(k) = 1;
                        end
                    end

                    if sum(f) == 0
                        check = 1;
                    end

                end

                stim_c = a(code_s);
                stimuli_design.(['design_' num2str(i)]).(['c_' num2str(cond(j))]) = stim_c;

            end

            if j == 3

                a = img_f.img_files_can_2;
                cond_pos = find(design_c(1,:)==cond(j));
                num_el = numel(cond_pos);
                rep = numel(cond_pos)/numel(a);

                code_c = 1:numel(a);
                code_r = repmat(code_c',rep,1);

                check = 0;
                while check == 0

                    code_s = Shuffle(code_r);
                    f = [];
                    for k = 1:length(code_c)
                        pat = findpattern(code_s,ones(2,1)*code_c(k));

                        if isempty(pat)
                            f(k) = 0;
                        else
                            f(k) = 1;
                        end
                    end

                    if sum(f) == 0
                        check = 1;
                    end

                end

                stim_c = a(code_s);
                stimuli_design.(['design_' num2str(i)]).(['c_' num2str(cond(j))]) = stim_c;

            end

            if j == 4

                a = img_f.img_files_noncan_2;
                cond_pos = find(design_c(1,:)==cond(j));
                num_el = numel(cond_pos);
                rep = numel(cond_pos)/numel(a);

                code_c = 1:numel(a);
                code_r = repmat(code_c',rep,1);

                check = 0;
                while check == 0

                    code_s = Shuffle(code_r);
                    f = [];
                    for k = 1:length(code_c)
                        pat = findpattern(code_s,ones(2,1)*code_c(k));

                        if isempty(pat)
                            f(k) = 0;
                        else
                            f(k) = 1;
                        end
                    end

                    if sum(f) == 0
                        check = 1;
                    end

                end

                stim_c = a(code_s);
                stimuli_design.(['design_' num2str(i)]).(['c_' num2str(cond(j))]) = stim_c;

            end

        end

        design_c = designs_exp.(['design_' num2str(i)]);
        stimuli_c = stimuli_design.(['design_' num2str(i)]);
        stimtype = design_c.design(1,:);

        o1_can_count = 1;
        o1_noncan_count = 1;
        o2_can_count = 1;
        o2_noncan_count = 1;

        stim_seq = cell(1,size(design_c.design,2));

        % Create stimulus vectors
        for j = 1:numel(stim_seq)

            if stimtype(j) == 1
                stim_seq(j) = stimuli_c.c_1(o1_can_count);
                o1_can_count = o1_can_count + 1;
            end
            if stimtype(j) == 2
                stim_seq(j) = stimuli_c.c_2(o1_noncan_count);
                o1_noncan_count = o1_noncan_count + 1;
            end
            if stimtype(j) == 3
                stim_seq(j) = stimuli_c.c_3(o2_can_count);
                o2_can_count = o2_can_count + 1;
            end
            if stimtype(j) == 4
                stim_seq(j) = stimuli_c.c_4(o2_noncan_count);
                o2_noncan_count = o2_noncan_count + 1;
            end
            if stimtype(j) == 0
                stim_seq(j) = {'NULL'};
            end

        end

        design_c.stim_seq = stim_seq;
        designs_exp.(['design_' num2str(i)]) = design_c;

        % Check design for repetitions between conditions
        stim_seq_c = designs_exp.(['design_' num2str(i)]).stim_seq;
        stim_seq_c = stim_seq_c';

        % Get rid of NULL trials
        idx = find(not(cellfun('isempty', strfind(stim_seq_c,'NULL'))));    
        stim_seq_c(idx) = [];

        stim_sequ_n = cell(length(stim_seq_c),1);

        for j = 1:length(stim_seq_c)

            name_c = stim_seq_c{j};
            us = strfind(name_c,'_');
            us = us(1);
            stim_sequ_n{j} = name_c(1:us-1);

        end

        stim_object = unique(stim_sequ_n);
        f = [];
        for j = 1:length(stim_object)

            stim_object_c = stim_object(j);
            stim_object_s = [stim_object_c;stim_object_c];

            d = [];
            for k = 1:length(stim_sequ_n)-1
                seq_stim = stim_sequ_n(k:k+1);

                diff = setdiff(seq_stim,stim_object_s);

                if isempty(diff)
                    d(k) = 1;
                else
                    d(k) = 0;
                end

            end

            f(j) = sum(d);

        end
        
        if sum(f) == 0
            check_r = 1;
            fprintf('Design %d done! \n',i)
            fprintf('Iterations: %d \n',count)
            fprintf('--------------- \n')
        end

    end
    
end

fprintf('Done! \n')
toc

save('Object_Can_NonCan_Localizer_Design_Stim.mat','designs_exp');
