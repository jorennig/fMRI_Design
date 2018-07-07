% Create designs from optseq par files
clc
clear all

% number of stimulus types
numstimtype = 3;
% stimulus names for filename output
stimtypes = {'MS1_F06_C1','MS1_F06_C3','MS1_F06_M','MS1_M03_C1','MS1_M03_C3','MS1_M03_M'};

data_files = dir('*.par');
for i=1:length(data_files)    
    design_files{i} = data_files(i).name;       
end
clear data_files

% number of runs
numruns = numel(design_files);
% number of stimuli per run
numstim = 90;
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
%out = 'McGurk_fMRI_ET';
designs = cell(numel(design_files),90);

for i = 1:numel(design_files)
    
    [stimonset,stimtype,stimlength,~,stimname] =  textread(design_files{i},'%f%n%f%f%s','delimiter',' ');
    designFile = [num2cell(stimonset),num2cell(stimtype),num2cell(stimlength),stimname];
    
    % We will have 2 seconds of data recording in the beginning
    stimonset = stimonset + 2;
    
    c1_1_onset = zeros(numel(find(stimtype==1)),1);
    c1_2_onset = zeros(numel(find(stimtype==2)),1);
    m1_onset = zeros(numel(find(stimtype==3)),1);
    c2_1_onset = zeros(numel(find(stimtype==4)),1);
    c2_2_onset = zeros(numel(find(stimtype==5)),1);
    m2_onset = zeros(numel(find(stimtype==6)),1);

    null_onset = zeros(numel(find(stimtype==0)),1);
    
    c1_1_duration = zeros(numel(find(stimtype==1)),1);
    c1_2_duration = zeros(numel(find(stimtype==2)),1);
    m1_duration = zeros(numel(find(stimtype==3)),1);
    c2_1_duration = zeros(numel(find(stimtype==4)),1);
    c2_2_duration = zeros(numel(find(stimtype==5)),1);
    m2_duration = zeros(numel(find(stimtype==6)),1);

    null_duration = zeros(numel(find(stimtype==0)),1);
    
    c1_1_count = 1;
    c1_2_count = 1;
    m1_count = 1;

    c2_1_count = 1;
    c2_2_count = 1;
    m2_count = 1;
    
    null_count = 1;
    
    % Create onset vectors
    for j = 1:numel(stimonset)

        if stimtype(j) == 1
            c1_1_onset(c1_1_count) = stimonset(j);
            c1_1_duration(c1_1_count) = stimlength(j);
            c1_1_count = c1_1_count + 1;
        end
        if stimtype(j) == 2
            c1_2_onset(c1_2_count) = stimonset(j);
            c1_2_duration(c1_2_count) = stimlength(j);
            c1_2_count = c1_2_count + 1;
        end
        if stimtype(j) == 3
            m1_onset(m1_count) = stimonset(j);
            m1_duration(m1_count) = stimlength(j);
            m1_count = m1_count + 1;
        end
        if stimtype(j) == 4
            c2_1_onset(c2_1_count) = stimonset(j);
            c2_1_duration(c2_1_count) = stimlength(j);
            c2_1_count = c2_1_count + 1;
        end
        if stimtype(j) == 5
            c2_2_onset(c2_2_count) = stimonset(j);
            c2_2_duration(c2_2_count) = stimlength(j);
            c2_2_count = c2_2_count + 1;
        end
        if stimtype(j) == 6
            m2_onset(m2_count) = stimonset(j);
            m2_duration(m2_count) = stimlength(j);
            m2_count = m2_count + 1;
        end        
        if stimtype(j) == 0
            null_onset(null_count) = stimonset(j);
            null_duration(null_count) = stimlength(j);
            null_count = null_count +1;
        end
                
    end
    
    % Save onsets of all conditions separately
    designs_exp.(['design_' num2str(i)]).onsets_mcgurk1 = m1_onset';
    designs_exp.(['design_' num2str(i)]).onsets_con1_1 = c1_1_onset';
    designs_exp.(['design_' num2str(i)]).onsets_con1_2 = c1_2_onset';
    designs_exp.(['design_' num2str(i)]).onsets_mcgurk2 = m2_onset';
    designs_exp.(['design_' num2str(i)]).onsets_con2_1 = c2_1_onset';
    designs_exp.(['design_' num2str(i)]).onsets_con2_2 = c2_2_onset';
    designs_exp.(['design_' num2str(i)]).onsets_null = null_onset';
    
    designs_exp.(['design_' num2str(i)]).duration_mcgurk1 = m1_duration';
    designs_exp.(['design_' num2str(i)]).duration_con1_1 = c1_1_duration';
    designs_exp.(['design_' num2str(i)]).duration_con1_2 = c1_2_duration';
    designs_exp.(['design_' num2str(i)]).duration_mcgurk2 = m2_duration';
    designs_exp.(['design_' num2str(i)]).duration_con2_1 = c2_1_duration';
    designs_exp.(['design_' num2str(i)]).duration_con2_2 = c2_2_duration';
    designs_exp.(['design_' num2str(i)]).duration_null = null_duration';

    % Save onsets for McGurk and Control trials
    m_onset = sort([m1_onset; m2_onset]);
    c1_onset = sort([c1_1_onset; c1_2_onset]);
    c2_onset = sort([c2_1_onset; c2_2_onset]);
    ba_onset = sort([c1_1_onset; c2_1_onset]);
    ga_onset = sort([c1_2_onset; c2_2_onset]);
    ma_onset = sort([m1_onset; c1_1_onset; c1_2_onset]);
    fe_onset = sort([m2_onset; c2_1_onset; c2_2_onset]);
    c_onset = sort([c1_1_onset; c1_2_onset; c2_1_onset; c2_2_onset]);

    designs_exp.(['design_' num2str(i)]).onsets_mcgurk = m_onset';
    designs_exp.(['design_' num2str(i)]).onsets_con1 = c1_onset';
    designs_exp.(['design_' num2str(i)]).onsets_con2 = c2_onset';
    designs_exp.(['design_' num2str(i)]).onsets_ba = ba_onset';
    designs_exp.(['design_' num2str(i)]).onsets_ga = ga_onset';
    
    designs_exp.(['design_' num2str(i)]).onsets_ma = ma_onset';
    designs_exp.(['design_' num2str(i)]).onsets_fe = fe_onset';

    designs_exp.(['design_' num2str(i)]).onsets_con = c_onset';
    
    designs_exp.(['design_' num2str(i)]).onsets_exp = stimonset';
    designs_exp.(['design_' num2str(i)]).stimduration_exp = stimlength';
    
    % optseq has the annoying habit of not listing multiple nulls in a row;
    % convert so each run has exactly the same number of stimuli
    design_exp = cell(numel(stimonset),1);

    for j = 1:numel(stimonset)

        if stimtype(j) == 0
            stim = 'NULL';

        elseif stimtype(j) == 1            
            stim = [stimtypes{1} '.mp4'];

        elseif stimtype(j) == 2
            stim = [stimtypes{2} '.mp4'];

        elseif stimtype(j) == 3
            stim = [stimtypes{3} '.mp4'];

        elseif stimtype(j) == 4            
            stim = [stimtypes{4} '.mp4'];

        elseif stimtype(j) == 5
            stim = [stimtypes{5} '.mp4'];

        elseif stimtype(j) == 6
            stim = [stimtypes{6} '.mp4'];

        end

        design_exp{j} = stim;

    end
    
    designs(i,1:numel(stimonset)) = design_exp;
    designs_exp.(['design_' num2str(i)]).design = [design_exp'; num2cell(stimonset)'];    
end % end reading in stimulus sequence fil

designs_exp.designs_overview = designs;

save('McGurk_fMRI_ET_2.mat','designs_exp');

% Create design txt files
for i = 1:numel(design_files)
    d = ['design_' num2str(i)];
    
    mc1(i,:) = designs_exp.(d).onsets_mcgurk1;
    c11(i,:) = designs_exp.(d).onsets_con1_1;
    c12(i,:) = designs_exp.(d).onsets_con1_2;

    mc2(i,:) = designs_exp.(d).onsets_mcgurk2;
    c21(i,:) = designs_exp.(d).onsets_con2_1;
    c22(i,:) = designs_exp.(d).onsets_con2_2;
    
    mc(i,:) = designs_exp.(d).onsets_mcgurk;
    c1(i,:) = designs_exp.(d).onsets_con1;
    c2(i,:) = designs_exp.(d).onsets_con2;
    ba(i,:) = designs_exp.(d).onsets_ba;
    ga(i,:) = designs_exp.(d).onsets_ga;
    
    ma(i,:) = designs_exp.(d).onsets_ma;
    fe(i,:) = designs_exp.(d).onsets_fe;
    
    c1and2(i,:) = designs_exp.(d).onsets_con;    
    
end

dlmwrite('MS1_F06_M.txt',mc1,'delimiter',' ','precision','%.1f');
dlmwrite('MS1_F06_C1.txt',c11,'delimiter',' ','precision','%.1f');
dlmwrite('MS1_F06_C3.txt',c12,'delimiter',' ','precision','%.1f');

dlmwrite('MS1_M03_M.txt',mc2,'delimiter',' ','precision','%.1f');
dlmwrite('MS1_M03_C1.txt',c21,'delimiter',' ','precision','%.1f');
dlmwrite('MS1_M03_C3.txt',c22,'delimiter',' ','precision','%.1f');

dlmwrite('MS1_M.txt',mc,'delimiter',' ','precision','%.1f');
dlmwrite('MS1_C1.txt',c1,'delimiter',' ','precision','%.1f');
dlmwrite('MS1_C3.txt',c2,'delimiter',' ','precision','%.1f');

dlmwrite('MS1_BA.txt',ba,'delimiter',' ','precision','%.1f');
dlmwrite('MS1_GA.txt',ga,'delimiter',' ','precision','%.1f');

dlmwrite('MS1_MA.txt',ma,'delimiter',' ','precision','%.1f');
dlmwrite('MS1_FE.txt',fe,'delimiter',' ','precision','%.1f');

dlmwrite('MS1_C.txt',c1and2,'delimiter',' ','precision','%.1f');

copyfile('MS1_F06_M.txt','/Users/beauchamplab/Documents/MATLAB/McGurk_Regressor_2');
copyfile('MS1_F06_C1.txt','/Users/beauchamplab/Documents/MATLAB/McGurk_Regressor_2');
copyfile('MS1_F06_C3.txt','/Users/beauchamplab/Documents/MATLAB/McGurk_Regressor_2');

copyfile('MS1_M03_M.txt','/Users/beauchamplab/Documents/MATLAB/McGurk_Regressor_2');
copyfile('MS1_M03_C1.txt','/Users/beauchamplab/Documents/MATLAB/McGurk_Regressor_2');
copyfile('MS1_M03_C3.txt','/Users/beauchamplab/Documents/MATLAB/McGurk_Regressor_2');

copyfile('MS1_M.txt','/Users/beauchamplab/Documents/MATLAB/McGurk_Regressor_2');
copyfile('MS1_C1.txt','/Users/beauchamplab/Documents/MATLAB/McGurk_Regressor_2');
copyfile('MS1_C3.txt','/Users/beauchamplab/Documents/MATLAB/McGurk_Regressor_2');

copyfile('MS1_BA.txt','/Users/beauchamplab/Documents/MATLAB/McGurk_Regressor_2');
copyfile('MS1_GA.txt','/Users/beauchamplab/Documents/MATLAB/McGurk_Regressor_2');

copyfile('MS1_MA.txt','/Users/beauchamplab/Documents/MATLAB/McGurk_Regressor_2');
copyfile('MS1_FE.txt','/Users/beauchamplab/Documents/MATLAB/McGurk_Regressor_2');

copyfile('MS1_C.txt','/Users/beauchamplab/Documents/MATLAB/McGurk_Regressor_2');
