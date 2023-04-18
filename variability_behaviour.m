function variability_behaviour(fig1)
%% variability_behaviour shows the variability in whisker behaviour
% This function plots:
% 1. The touch epochs and delta kappa from an example trial
% 2. The touch epochs from an example session
% 3. Statistics of behavioural variables (divided into hit and miss trials):
%   3.1 Delta Kappa 
%   3.2 Touch length
%   3.3 Inter-touch interval
%   3.4 Number of touches

ff=dir('*.mat*');
% LT=[];
% NT=[];
% ITT=[];
% dk=[];
dk_correct=[];
dk_incorrect=[];
LT_correct=[];
LT_incorrect=[];
NT_correct=[];
NT_incorrect=[];
ITT_correct=[];
ITT_incorrect=[];
percentile_1=1;
percentile_2=1;
figure(fig1)
for f=1:size(ff,1)
    load(ff(f).name,'Data')
    Ntrials=size(Data.touch,1);
    total_psth=ones(size(Data.touch));


    if strcmp(ff(f).name,'Glu43_22122017H.mat')

        subplot(4,2,[ 3 5 7])
        [~]=touch_variability_behaviour(Data,1:Ntrials,1);
 
        %[nt_per_trial,inter_time_touch,length_touch]=touch_variability_behaviour(Data,1:Ntrials,1);
    end
%     else
%         [nt_per_trial,inter_time_touch,length_touch]=touch_variability_behaviour(Data,1:Ntrials,0);
%     end
    
    %[~,~,~,~,~,~,deltak]=Tcurve_touches_dk(Data.touch,Data.touch_per_whisker,total_psth,Data.deltak_w1,Data.deltak_w2,0,4,0,percentile_1,percentile_2);
    
    % only correct trials
    [~,~,~,~,~,~,deltak_correct]=Tcurve_touches_dk(Data.touch(Data.correct_trials,:,:),Data.touch_per_whisker(Data.correct_trials,:,:),total_psth(Data.correct_trials,:),Data.deltak_w1(Data.correct_trials,:),Data.deltak_w2(Data.correct_trials,:),0,4,0,percentile_1,percentile_2);
    [nt_per_trial_correct,inter_time_touch_correct,length_touch_correct]=touch_variability_behaviour(Data,Data.correct_trials,0);
    
    % only incorrect trials
    [~,~,~,~,~,~,deltak_incorrect]=Tcurve_touches_dk(Data.touch(Data.incorrect_trials,:,:),Data.touch_per_whisker(Data.incorrect_trials,:,:),total_psth(Data.incorrect_trials,:,:),Data.deltak_w1(Data.incorrect_trials,:,:),Data.deltak_w2(Data.incorrect_trials,:,:),0,4,0,percentile_1,percentile_2);
    [nt_per_trial_incorrect,inter_time_touch_incorrect,length_touch_incorrect]=touch_variability_behaviour(Data,Data.incorrect_trials,0);
    
    
    
%    dk=[dk;deltak'];
%     LT=[LT;length_touch'];
%     NT=[NT;nt_per_trial];
%     ITT=[ITT;inter_time_touch'];


    dk_correct=[dk_correct;deltak_correct'];
    dk_incorrect=[dk_incorrect;deltak_incorrect'];
    LT_correct=[LT_correct;length_touch_correct'];
    LT_incorrect=[LT_incorrect;length_touch_incorrect'];
    NT_correct=[NT_correct;nt_per_trial_correct];
    NT_incorrect=[NT_incorrect;nt_per_trial_incorrect];
    ITT_correct=[ITT_correct;inter_time_touch_correct'];
    ITT_incorrect=[ITT_incorrect;inter_time_touch_incorrect'];
    
    clear inter_time_touch_correct length_touch_correct length_touch_incorrect nt_per_trial_correct nt_per_trial_incorrect length_touch nt_per_trial inter_time_touch deltak go_trials info idx_correct idx_incorrect deltak_correct deltak_incorrect
    
end

subplot(4,2,2)
hold on
histogram(dk_correct,-0.03:0.001:0.03,'Normalization','Probability','FaceColor','g')
histogram(dk_incorrect,-0.03:0.001:0.03,'Normalization','Probability','FaceColor', 'r')
box off
xticks([-0.03 -0.02 -0.01 0 0.01 0.02 0.03])
yticks([0 0.1 0.2 0.3])
ylim([0 0.3])
xlabel('\Delta \kappa [mm^{-1}]')
ylabel('Fraction of touches')

subplot(4,2,4)
histogram(LT_correct,0:10:350,'Normalization','Probability','FaceColor','g')
hold on
histogram(LT_incorrect,0:10:350,'Normalization','Probability','FaceColor','r')
xticks([0 100 200 300])
yticks([0 0.1 0.2 0.3])
box off
xlabel('Touch  length  [ms]')
ylabel('Fraction of touches')


% inset

ax1=axes('Position',[0.8 0.6 0.1 0.1]);
histogram(ax1,LT_correct,1000:100:3000,'Normalization','Probability','FaceColor','g')
hold on
histogram(ax1,LT_incorrect,1000:100:3000,'Normalization','Probability','FaceColor','r')
box off
xticks([1000 3000])
yticks([0 0.05])




subplot(4,2,6)
histogram(ITT_correct,0:10:300,'Normalization','Probability','FaceColor','g')
hold on
histogram(ITT_incorrect,0:10:300,'Normalization','Probability','FaceColor','r')
box off
xticks([0 100 200 300])
yticks([0 0.1 0.2 0.3])
ylim([0 0.3])
xlabel('Inter-touch interval [ms]')
ylabel('Fraction of touches')

subplot(4,2,8)
box off
histogram(NT_correct,0:1:40,'Normalization','Probability','FaceColor','g')
hold on
histogram(NT_incorrect,0:1:40,'Normalization','Probability','FaceColor','r')
xticks([0 10 20 30 40])
yticks([0 0.1 0.2 0.3])
ylim([0 0.3])
xlabel('Number of touches')
ylabel('Fraction of touches')
legend('Hit trials','Miss trials')
box off

%% Averages
% average number of touches in hit trials
disp(['Mean number of touches hit trials ' num2str(mean(NT_correct)) '+-' num2str(std(NT_correct))])
%average number of touches in miss trials
disp(['Mean number of touches miss trials ' num2str(mean(NT_incorrect)) '+-' num2str(std(NT_incorrect))])
end