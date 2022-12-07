function variability_behaviour
folder=[{'Glu32_19092017H'};{'Glu32_21092017H'};{'Glu43_22122017H'};{'Glu35_10112017H'};{'Glu35_13112017H_1'};{'Glu35_13112017H_2'}];
LT=[];
NT=[];
ITT=[];
dk=[];
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
for f=1:size(folder,1)
    cd([folder{f} '/'])
    info=xlsread('ledtrials.xlsx');
    go_trials=find(info(:,7)==2);
    idx_correct=info(go_trials,6)==1;
    idx_incorrect=info(go_trials,6)==2;
    
    load('touches_whisker.mat','touches_whisker')
    total_psth=ones(size(touches_whisker(:,:,1)));
    load('newk1.mat','newk1')
    k_c1=newk1*1.7/0.047; %%corrected by thickness
    load('newk2.mat','newk2')
    k_c2=newk2/0.047;
    load('all_touches.mat','touches_matrix')
    
    % original line (all trials)
    if f==3
        [nt_per_trial,inter_time_touch,length_touch]=touch_variability_behaviour(1:size(k_c2,1),1);
    else
        [nt_per_trial,inter_time_touch,length_touch]=touch_variability_behaviour(1:size(k_c2,1),0);
    end
    
    [~,~,~,~,~,~,~,~,~,~,deltak,~,~,~,~,~]=Tcurve_touches_dk(touches_matrix,touches_whisker,total_psth,k_c1,k_c2,0,4,0,percentile_1,percentile_2);
    % only correct trials
    [~,~,~,~,~,~,~,~,~,~,deltak_correct,~,~,~,~,~]=Tcurve_touches_dk(touches_matrix(idx_correct,:,:),touches_whisker(idx_correct,:,:),total_psth(idx_correct,:,:),k_c1(idx_correct,:,:),k_c2(idx_correct,:,:),0,4,0,percentile_1,percentile_2);
    [nt_per_trial_correct,inter_time_touch_correct,length_touch_correct]=touch_variability_behaviour(idx_correct,0);
    % only incorrect trials
    [~,~,~,~,~,~,~,~,~,~,deltak_incorrect,~,~,~,~,~]=Tcurve_touches_dk(touches_matrix(idx_incorrect,:,:),touches_whisker(idx_incorrect,:,:),total_psth(idx_incorrect,:,:),k_c1(idx_incorrect,:,:),k_c2(idx_incorrect,:,:),0,4,0,percentile_1,percentile_2);
    [nt_per_trial_incorrect,inter_time_touch_incorrect,length_touch_incorrect]=touch_variability_behaviour(idx_incorrect,0);
    
    
    
    dk=[dk;deltak'];
    dk_correct=[dk_correct;deltak_correct'];
    dk_incorrect=[dk_incorrect;deltak_incorrect'];
    LT=[LT;length_touch'];
    NT=[NT;nt_per_trial];
    ITT=[ITT;inter_time_touch'];
    
    LT_correct=[LT_correct;length_touch_correct'];
    LT_incorrect=[LT_incorrect;length_touch_incorrect'];
    NT_correct=[NT_correct;nt_per_trial_correct];
    NT_incorrect=[NT_incorrect;nt_per_trial_incorrect];
    ITT_correct=[ITT_correct;inter_time_touch_correct'];
    ITT_incorrect=[ITT_incorrect;inter_time_touch_incorrect'];
    
    clear inter_time_touch_correct length_touch_correct length_touch_incorrect nt_per_trial_correct nt_per_trial_incorrect length_touch nt_per_trial inter_time_touch deltak go_trials info idx_correct idx_incorrect deltak_correct deltak_incorrect
    cd ..
    
end
figure
subplot(4,1,1)
hold on
histogram(dk_correct,-0.03:0.001:0.03,'Normalization','Probability','FaceColor','g')
histogram(dk_incorrect,-0.03:0.001:0.03,'Normalization','Probability','FaceColor', 'r')
box off
xticks([-0.03 -0.02 -0.01 0 0.01 0.02 0.03])
yticks([0 0.1 0.2 0.3])
ylim([0 0.3])
xlabel('\Delta \kappa [mm^{-1}]')
ylabel('Fraction of touches')

subplot(4,1,2)
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
histogram(LT_correct,1000:100:3000,'Normalization','Probability','FaceColor','g')
hold on
histogram(LT_incorrect,1000:100:3000,'Normalization','Probability','FaceColor','r')
box off
xticks([1000 3000])
yticks([0 0.05])




subplot(4,1,3)
histogram(ITT_correct,0:10:300,'Normalization','Probability','FaceColor','g')
hold on
histogram(ITT_incorrect,0:10:300,'Normalization','Probability','FaceColor','r')
box off
xticks([0 100 200 300])
yticks([0 0.1 0.2 0.3])
ylim([0 0.3])
xlabel('Inter-touch interval [ms]')
ylabel('Fraction of touches')

subplot(4,1,4)
box off
histogram(NT_correct,0:1:40,'Normalization','Probability','FaceColor','g')
hold on
histogram(NT_incorrect,0:1:40,'Normalization','Probability','FaceColor','r')


%% averages
% average number of touches in hit trials
mean_std_Number_of_T_hit_trials=[mean(NT_correct) std(NT_correct)]
%average number of touches in miss trials
mean_std_Number_of_T_miss_trials=[mean(NT_incorrect) std(NT_incorrect)]
xticks([0 10 20 30 40])
yticks([0 0.1 0.2 0.3])
ylim([0 0.3])
xlabel('Number of touches')
ylabel('Fraction of touches')
legend('Hit trials','Miss trials')
box off
end