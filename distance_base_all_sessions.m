function distance_base_all_sessions(fig3)
ff=dir('*.mat*');
nfolders=size(ff,1);
correct=[];
y=[];
first_later=[];
per_session=nan(nfolders,4);
per_session_correct=nan(nfolders,4);
per_session_incorrect=nan(nfolders,4);

figure(fig3)
subplot(4,3,8)
hold on 

for f=1:nfolders
    
    load(ff(f).name,'Data')
    
    [distance_1,distance_2,distance_3,distance_4,y_tmp,touch_idx,hit_tmp]=distance_base_from_matrix(Data,ones(size(Data.correct_trials))>0);
    [d1_correct,d2_correct,d3_correct,d4_correct,~,~,~]=distance_base_from_matrix(Data,Data.correct_trials);
    [d1_incorrect,d2_incorrect,d3_incorrect,d4_incorrect,~,~,~]=distance_base_from_matrix(Data,Data.incorrect_trials);
    
    per_session(f,:)=[mean(distance_1) mean(distance_2) mean(distance_3) mean(distance_4)]*0.047;
    
    
    per_session_correct(f,:)=[mean(d1_correct) mean(d2_correct) mean(d3_correct) mean(d4_correct)]*0.047;
    
    
    per_session_incorrect(f,:)=[mean(d1_incorrect) mean(d2_incorrect) mean(d3_incorrect) mean(d4_incorrect)]*0.047;
    
    % Session 3 is clearly an outlier, so it's removed from the statistical
    % test
    
    if ~strcmp('Glu43_22122017H.mat',ff(f).name)
        y=[y;y_tmp(touch_idx==1)';y_tmp(touch_idx~=1)'];
        first_later=[first_later;touch_idx(touch_idx==1)';touch_idx(touch_idx~=1)'*0];
        correct=[correct;hit_tmp(touch_idx==1)';hit_tmp(touch_idx~=1)'];
    end
    
    %[sum(hit_tmp'==1) sum(hit_tmp'==0 & touch_idx'==1)  sum(touch_idx'==1) sum(touch_idx'~=1)]
    %p = anovan(y_tmp',{hit_tmp',touch_idx'==1});
    
 
    
    plot(per_session(f,:),'Color',[0.5 0.5 0.5])
    
end


hold on
errorbar(mean(per_session,'omitnan'),std(per_session,'omitnan'),'k.-')
errorbar(mean(per_session_correct,'omitnan'),std(per_session_correct,'omitnan'),'g.-')
errorbar(mean(per_session_incorrect,'omitnan'),std(per_session_incorrect,'omitnan'),'r.-')

xlabel('N touches')
ylabel('Speed [mm/ms]')
ylim([0 0.04])
box off


%ANOVA test for figure 3
y=y*0.047; %convert pixels to mm
% unbalanced
p = anovan(y,{correct,first_later},'display','off');

disp(' ')
disp('-----------------------------------')
disp('Unbalanced anova First vs later touches speed')
disp(['p-value Effect hit vs miss = ' num2str(p(1))])
disp(['p-value Effect first vs later = ' num2str(p(2))])

% balanced
nsample=sum(correct==0 & first_later==1);
correct_1=y(correct==1 & first_later==1);
incorrect_1=y(correct==0 & first_later==1);
correct_later=y(correct==1 & first_later==0);
incorrect_later=y(correct==0 & first_later==0);
balanced_anova(correct_1,incorrect_1,correct_later,incorrect_later,nsample)

% proportion_correct_incorrect=mean(y(correct==1))/mean(y(correct==0))
% proportion_first_later=mean(y(first_later~=1))/mean(y(first_later==1))

end