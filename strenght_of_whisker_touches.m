function strenght_of_whisker_touches(fig3)
folder=[{'Glu32_19092017H'};{'Glu32_21092017H'};{'Glu43_22122017H'};{'Glu35_10112017H'};{'Glu35_13112017H_1'};{'Glu35_13112017H_1'}];
nfolders=size(folder,1);
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
for f=1:nfolders
    
    cd(folder{f})
    info=xlsread('ledtrials.xlsx');
    go_trials=find(info(:,7)==2);
    idx_correct=info(go_trials,6)==1;
    idx_incorrect=info(go_trials,6)==2;
    
    
    load('newk1.mat','newk1')
    k_c1=newk1*1.7/0.047;
    load('newk2.mat','newk2')
    k_c2=newk2/0.047;
    
    load('all_touches.mat','touches_matrix')
    load('touches_whisker.mat','touches_whisker')
    
    %% select correct and incorrect trials
    touches_matrix_correct=touches_matrix(idx_correct,:);
    touches_whisker_correct=touches_whisker(idx_correct,:,:);
    k_c1_correct=k_c1(idx_correct,:);
    k_c2_correct=k_c2(idx_correct,:);
    touches_matrix_incorrect=touches_matrix(idx_incorrect,:);
    touches_whisker_incorrect=touches_whisker(idx_incorrect,:,:);
    k_c1_incorrect=k_c1(idx_incorrect,:);
    k_c2_incorrect=k_c2(idx_incorrect,:);
    
    total_psth=touches_matrix*0+10;
    total_psth_correct=touches_matrix*0+10;
    total_psth_incorrect=touches_matrix*0+10;
    
    [~,~,~,~,~,~,~,~,~,~,deltak,touch_idx,~,~,~,~]=Tcurve_touches_dk(touches_matrix,touches_whisker,total_psth,k_c1,k_c2,0,4,0,percentile_1,percentile_2);
    [~,~,~,~,~,~,~,~,~,~,deltak_correct,touch_idx_correct,~,~,~,~]=Tcurve_touches_dk(touches_matrix_correct,touches_whisker_correct,total_psth_correct,k_c1_correct,k_c2_correct,0,4,0,percentile_1,percentile_2);
    [~,~,~,~,~,~,~,~,~,~,deltak_incorrect,touch_idx_incorrect,~,~,~,~]=Tcurve_touches_dk(touches_matrix_incorrect,touches_whisker_incorrect,total_psth_incorrect,k_c1_incorrect,k_c2_incorrect,0,4,0,percentile_1,percentile_2);
    
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
    cd ..
    
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

[hd1,pd1] = ttest2(deltak_touch(:,1),deltak_touch(:,2));
[hd2,pd2] = ttest2(deltak_touch(:,2),deltak_touch(:,3));
[hd3,pd3] = ttest2(deltak_touch(:,3),deltak_touch(:,4));



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


%mixing all trials
[h_dk_1,pdk_1] = ttest2(deltak_correct_1_all,deltak_incorrect_1_all);
[h_dk_later,pdk_later] = ttest2(deltak_correct_later_all,deltak_incorrect_later_all);

[h_dk_correct,pdk_correct] = ttest2(deltak_correct_1_all,deltak_correct_later_all);
[h_dk_incorrect,pdk_incorrect] = ttest2(deltak_incorrect_1_all,deltak_incorrect_later_all);

%%%% ANOVA figure 3
%% unbalanced (raw)

correct=[ones(size(deltak_correct_1_all));zeros(size(deltak_incorrect_1_all));ones(size(deltak_correct_later_all));zeros(size(deltak_incorrect_later_all))];
first_later=[ones(size(deltak_correct_1_all));ones(size(deltak_incorrect_1_all));ones(size(deltak_correct_later_all))*2;ones(size(deltak_incorrect_later_all))*2];
y=[deltak_correct_1_all;deltak_incorrect_1_all;deltak_correct_later_all;deltak_incorrect_later_all];
p = anovan(y,{correct,first_later});
p(1)
p(2)
%% balanced
%%minimun number of samples per group: size(deltak_incorrect_1_all)
% for i=1:100
% nsample=size(deltak_incorrect_1_all,1);
% correct=[ones(nsample,1);zeros(nsample,1);ones(nsample,1);zeros(nsample,1)];
% first_later=[ones(nsample,1);ones(nsample,1);ones(nsample,1)*2;ones(nsample,1)*2];
% y=[deltak_correct_1_all(randperm(numel(deltak_correct_1_all),nsample));...
%     deltak_incorrect_1_all(randperm(numel(deltak_incorrect_1_all),nsample));...
%     deltak_correct_later_all((randperm(numel(deltak_correct_later_all),nsample)));...
%     deltak_incorrect_later_all((randperm(numel(deltak_incorrect_later_all),nsample)))];
% p(:,i) = anovan(y,{correct,first_later},'display','off');
%
% end

end