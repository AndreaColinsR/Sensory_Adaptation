function strength_of_whisker_touches(fig3)
%% strength_of_whisker_touches plots 
% 1. Delta kappa per touch number divided into hit and miss trials 
% 2. Delta kappa per session divided into first/ later and hit/miss trials

% this function also does ANOVA (balanced and unbalanced) for delta kappa for first/ later and
% hit/miss trials conditions

ff=dir('*.mat*');
nfolders=size(ff,1);
percentile_1=0.95;
percentile_2=0.95;

deltak_touch=nan(nfolders,4);
deltak_touch_correct=nan(nfolders,4);
deltak_touch_incorrect=nan(nfolders,4);
deltak_correct_incorrect_1=nan(nfolders,2);
deltak_correct_incorrect_later=nan(nfolders,2);

deltak_correct_1_all=[];
deltak_incorrect_1_all=[];
deltak_correct_later_all=[];
deltak_incorrect_later_all=[];

figure(fig3)

for f=1:size(ff,1)
    

    load(ff(f).name,'Data')
    
    
    %% select correct and incorrect trials
    touches_matrix_correct=Data.touch(Data.correct_trials,:);
    touches_whisker_correct=Data.touch_per_whisker(Data.correct_trials,:,:);
    k_c1_correct=Data.deltak_w1(Data.correct_trials,:);
    k_c2_correct=Data.deltak_w2(Data.correct_trials,:);
    
    touches_matrix_incorrect=Data.touch(Data.incorrect_trials,:);
    touches_whisker_incorrect=Data.touch_per_whisker(Data.incorrect_trials,:,:);
    k_c1_incorrect=Data.deltak_w1(Data.incorrect_trials,:);
    k_c2_incorrect=Data.deltak_w2(Data.incorrect_trials,:);
    
    total_psth=Data.touch*0+10;
    total_psth_correct=Data.touch*0+10;
    total_psth_incorrect=Data.touch*0+10;
    
    [~,~,~,~,~,~,deltak,touch_idx]=Tcurve_touches_dk(Data.touch,Data.touch_per_whisker,total_psth,Data.deltak_w1,Data.deltak_w2,0,4,0,percentile_1,percentile_2);
    [~,~,~,~,~,~,deltak_correct,touch_idx_correct]=Tcurve_touches_dk(touches_matrix_correct,touches_whisker_correct,total_psth_correct,k_c1_correct,k_c2_correct,0,4,0,percentile_1,percentile_2);
    [~,~,~,~,~,~,deltak_incorrect,touch_idx_incorrect]=Tcurve_touches_dk(touches_matrix_incorrect,touches_whisker_incorrect,total_psth_incorrect,k_c1_incorrect,k_c2_incorrect,0,4,0,percentile_1,percentile_2);
    

    deltak=abs(deltak);
    deltak_correct=abs(deltak_correct);
    deltak_incorrect=abs(deltak_incorrect);

    
    deltak_touch(f,:)=[mean(deltak(touch_idx==1)) mean(deltak(touch_idx==2)) mean(deltak(touch_idx==3)) mean(deltak(touch_idx>3))];
    deltak_touch_correct(f,:)=[mean(deltak_correct(touch_idx_correct==1)) mean(deltak_correct(touch_idx_correct==2)) mean(deltak_correct(touch_idx_correct==3)) mean(deltak_correct(touch_idx_correct>3))];
    deltak_touch_incorrect(f,:)=[mean(deltak_incorrect(touch_idx_incorrect==1)) mean(deltak_incorrect(touch_idx_incorrect==2)) mean(deltak_incorrect(touch_idx_incorrect==3)) mean(deltak_incorrect(touch_idx_incorrect>3))];
    
    
    deltak_correct_incorrect_1(f,:)=[mean(deltak_correct(touch_idx_correct==1)) mean(deltak_incorrect(touch_idx_incorrect==1))];
    deltak_correct_incorrect_later(f,:)=[mean(deltak_correct(touch_idx_correct>1)) mean(deltak_incorrect(touch_idx_incorrect>1))];
    
    deltak_correct_1_all=[deltak_correct_1_all;deltak_correct(touch_idx_correct==1)'];
    deltak_incorrect_1_all=[deltak_incorrect_1_all;deltak_incorrect(touch_idx_incorrect==1)'];
    deltak_correct_later_all=[deltak_correct_later_all;deltak_correct(touch_idx_correct>1)'];
    deltak_incorrect_later_all=[deltak_incorrect_later_all;deltak_incorrect(touch_idx_incorrect>1)'];
     
    
    
    clear deltak touch_idx

end

subplot(4,3,7)
hold on
plot(1:4,deltak_touch,'Color',[0.5 0.5 0.5])
xlabel('Touch number')
xticks([1 2 3 4])
ylabel('\Delta \kappa_{3D} [1/mm]')
xlim([0.5 4.5])

errorbar(1:4,mean(deltak_touch,1),std(deltak_touch,[],1),'k')
errorbar(1:4,mean(deltak_touch_correct,1),std(deltak_touch_correct,[],1),'g')
errorbar(1:4,mean(deltak_touch_incorrect,1),std(deltak_touch_incorrect,[],1),'r')


subplot(4,3,9)
hold on
errorbar([0.8 1.8],[mean(deltak_correct_1_all) mean(deltak_incorrect_1_all)],[std(deltak_correct_1_all) std(deltak_incorrect_1_all)],'.m')
errorbar([1.2 2.2],mean(deltak_correct_incorrect_later),std(deltak_correct_incorrect_later),'.b')

plot([0.8 1.8],deltak_correct_incorrect_1,'.r')
plot([1.2 2.2],deltak_correct_incorrect_later,'.k')

ylabel('\Delta \kappa_{3D} [1/mm]')
xlim([0.6 2.4])
box off
xticks([1 2])
xticklabels({'Hit','Miss'})


%%%% ANOVA figure 3
%% unbalanced (raw)
disp('Delta kappa comparison')
correct=[ones(size(deltak_correct_1_all));zeros(size(deltak_incorrect_1_all));ones(size(deltak_correct_later_all));zeros(size(deltak_incorrect_later_all))];
first_later=[ones(size(deltak_correct_1_all));ones(size(deltak_incorrect_1_all));ones(size(deltak_correct_later_all))*2;ones(size(deltak_incorrect_later_all))*2];

y=[deltak_correct_1_all;deltak_incorrect_1_all;deltak_correct_later_all;deltak_incorrect_later_all];
p = anovan(y,{correct,first_later},'display','off');

disp('-----------------------------------')
disp('Unbalanced anova First vs later Delta Kappa')
disp(['p-value Effect hit vs miss = ' num2str(p(1))])
disp(['p-value Effect first vs later = ' num2str(p(2))])

%% balanced
nsample=size(deltak_incorrect_1_all,1);
balanced_anova(deltak_correct_1_all,deltak_incorrect_1_all,deltak_correct_later_all,deltak_incorrect_later_all,nsample)
end