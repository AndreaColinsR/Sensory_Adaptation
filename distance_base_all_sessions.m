function distance_base_all_sessions(fig3)
folder=[{'Glu32_19092017H'};{'Glu32_21092017H'};{'Glu43_22122017H'};{'Glu35_10112017H'};{'Glu35_13112017H_1'};{'Glu35_13112017H_2'}];
nfolders=size(folder,1);
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
    cd(folder{f})
    info=xlsread('ledtrials.xlsx');
    idx_all=info(:,7)==2;
    idx_correct=info(idx_all,6)==1;
    idx_incorrect=info(idx_all,6)==2;
    
    %[distance_1,distance_2,distance_3,distance_4,distance_matrix,y_tmp,touch_idx,hit_tmp]=distance_base(folder2{f},idx_all);
    %[d1_correct,d2_correct,d3_correct,d4_correct,~,~,~,~]=distance_base(folder2{f},idx_correct);
    %[d1_incorrect,d2_incorrect,d3_incorrect,d4_incorrect,~,~,~,~]=distance_base(folder2{f},idx_incorrect);
    
    [distance_1,distance_2,distance_3,distance_4,y_tmp,touch_idx,hit_tmp]=distance_base_from_matrix(idx_all);
    [d1_correct,d2_correct,d3_correct,d4_correct,~,~,~]=distance_base_from_matrix(idx_correct);
    [d1_incorrect,d2_incorrect,d3_incorrect,d4_incorrect,~,~,~]=distance_base_from_matrix(idx_incorrect);
    
    per_session(f,:)=[mean(distance_1) mean(distance_2) mean(distance_3) mean(distance_4)]*0.047;
    
    
    per_session_correct(f,:)=[mean(d1_correct) mean(d2_correct) mean(d3_correct) mean(d4_correct)]*0.047;
    
    
    per_session_incorrect(f,:)=[mean(d1_incorrect) mean(d2_incorrect) mean(d3_incorrect) mean(d4_incorrect)]*0.047;
    
    % Session 3 is clearly an outlier, so it's removed from the statistical
    % test
    
    if f~=3
        y=[y;y_tmp(touch_idx==1)';y_tmp(touch_idx~=1)'];
        first_later=[first_later;touch_idx(touch_idx==1)';touch_idx(touch_idx~=1)'*0];
        correct=[correct;hit_tmp(touch_idx==1)';hit_tmp(touch_idx~=1)'];
    end
    %[sum(hit_tmp'==1) sum(hit_tmp'==0 & touch_idx'==1)  sum(touch_idx'==1) sum(touch_idx'~=1)]
    %p = anovan(y_tmp',{hit_tmp',touch_idx'==1});
    
    
    %     y=[y;y_tmp'];
    %     first_later=[first_later;touch_idx'==1];
    %     correct=[correct;hit_tmp'];
    
    plot(per_session(f,:),'Color',[0.5 0.5 0.5])
    
    cd ..
end


hold on
errorbar(mean(per_session,'omitnan'),std(per_session,'omitnan'),'k.-')
errorbar(mean(per_session_correct,'omitnan'),std(per_session_correct,'omitnan'),'g.-')
errorbar(mean(per_session_incorrect,'omitnan'),std(per_session_incorrect,'omitnan'),'r.-')

xlabel('N touches')
ylabel('Speed [mm/ms]')
ylim([0 0.04])
box off

% t test
% [h1,p1] = ttest(per_session(:,1),per_session(:,2),'Tail','both')
% [h2,p2] = ttest(per_session(:,2),per_session(:,3),'Tail','both')
% [h3,p3] = ttest(per_session(:,3),per_session(:,4),'Tail','both')

%ANOVA test for figure 3
y=y*0.047; %convert pixels to mm
idx=y<1;

p = anovan(y(idx),{correct(idx),first_later(idx)});
p(1)
p(2)


% figure
% subplot(2,1,1)
% errorbar([1 2],[mean(y(correct==1)) mean(y(correct==0))],[std(y(correct==1)) std(y(correct==0))])


% hold on
% plot(1,y(correct==1),'.')
% plot(2,y(correct==0),'.')
%y_correct=y(correct==1);
%X = randi(numel(y_correct),[sum(correct==0) 1]);
% ttest(y_correct(X),y(correct==0))
% ttest(y(correct==1),y(correct==0))
% histogram(y_correct(X),0:0.01:0.3)
% hold on
% histogram(y(correct==0),0:0.01:0.3)

% subplot(2,1,2)
% errorbar([1 2],[mean(y(first_later==1)) mean(y(first_later~=1))],[std(y(first_later==1)) std(y(first_later~=1))])

proportion_correct_incorrect=mean(y(correct==1))/mean(y(correct==0))
proportion_first_later=mean(y(first_later~=1))/mean(y(first_later==1))
%[sum(correct==1) sum(correct==0) sum(first_later==1) sum(first_later~=1)]

% box off
% xticks([1 2 3 4])
% ylim([0 0.04])
% xlim([0.8 4.2])
% yticks([0 0.01 0.02 0.03 0.04])
end