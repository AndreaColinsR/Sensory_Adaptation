function change_format_all_datasets
all_folders=dir('Glu*');
%total_units=0;
for f=1:size(all_folders,1)
    f
    cd(all_folders(f).name)
    [deltak_w1,deltak_w2,touch,touch_per_whisker,correct_trials,incorrect_trials,unit]=Change_format_dataset;
    %% chheck that the total number of units is correct N=33;
    %total_units=total_units+size(unit,2);
    cd ..
    save([all_folders(f).name '.mat'],'deltak_w1','deltak_w2','touch','touch_per_whisker','correct_trials','incorrect_trials','unit')
    
    clear deltak_w1 deltak_w2 touch touch_per_whisker correct_trials incorrect_trials unit
end

end


function [deltak_w1,deltak_w2,touch,touch_per_whisker,correct_trials,incorrect_trials,unit]=Change_format_dataset
%% Load all the files in the folder and put them in the correct
%format

%% Behaviour


% this is Delta K
load('newk1.mat')
deltak_w1=newk1*1.7/0.047;
load('newk2.mat')
deltak_w2=newk2/0.047;

% this is touch independent of the whisker
load('all_touches.mat')
touch=touches_matrix;

% this is touch independent of the whisker
load('touches_whisker.mat')
touch_per_whisker=touches_whisker;


% correct and incorrect trials
trial_output=table2array(readtable('ledtrials.xlsx','Range','F:G'));


go_trial=find(trial_output(:,2)==2);

correct_trials=zeros(size(go_trial,1),1);
incorrect_trials=zeros(size(go_trial,1),1);

correct_trials(trial_output(go_trial,1)==1)=1;
incorrect_trials(trial_output(go_trial,1)==2)=1;

% check that this is the same than what we have before
info=xlsread('ledtrials.xlsx');
go_trials=find(info(:,7)==2);
idx_correct=info(go_trials,6)==1;
idx_incorrect=info(go_trials,6)==2;

if any(correct_trials~=idx_correct)
    disp('Correct trials are not the same')
    keyboard
end

if any(incorrect_trials~=idx_incorrect)
    disp('Correct trials are not the same')
    keyboard
end


%% Neural activity
neurons=dir('neural_data*');
for in=1:size(neurons,1)
    load(neurons(in).name,'psth')
    unit(in).spikes=psth;
    clear psth
end

end