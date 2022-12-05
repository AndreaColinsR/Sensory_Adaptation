function [AUC_1,AUC_later,AUC_12,AUC_later2]=recovery_test
figure
folder=[{'Glu32_19092017H'};{'Glu32_21092017H'};{'Glu43_22122017H'};{'Glu35_10112017H'};{'Glu35_13112017H/first_segment'};{'Glu35_13112017H/second_segment'}];
counter=1;
[~,txt]=xlsread('./Adaptation/Adaptation_units_list.xlsx',1);
for f=1:size(folder,1)
    cd(folder{f})
    info=xlsread('ledtrials.xlsx');
    go_trials=find(info(:,7)==2);
    idx_correct=info(go_trials,6)==1;
    idx_incorrect=info(go_trials,6)==2;
    cd('./modelling/')
    
    ff=dir('*neural_data_S*');
    load(ff(1).name)
    total_psth=zeros(size(psth));
    for unit=1:size(ff,1)
        
         if onthelist(txt,[folder{f} ff(unit).name])
        load(ff(unit).name)
        total_psth=psth;
   
    
    [AUC_1_correct(counter),AUC_later_correct(counter),AUC_12_correct(counter),AUC_later2_correct(counter)]=Detection_test(total_psth(idx_correct,:),idx_correct);

    [AUC_1_incorrect(counter),AUC_later_incorrect(counter),AUC_12_incorrect(counter),AUC_later2_incorrect(counter)]=Detection_test(total_psth(idx_incorrect,:),idx_incorrect);
    %pause
    [AUC_1(counter),AUC_later(counter),AUC_12(counter),AUC_later2(counter)]=Detection_test(total_psth);

    counter=counter+1;
    %pause
    clear total_psth
         end
    end
    if f==size(folder,1)-1
        cd ..
        cd ..
        cd ..
    else
        cd ..
        cd ..
        
    end
end
subplot(2,2,4)
errorbar(1:4,[mean(AUC_1_correct) mean(AUC_later_correct) mean(AUC_1_incorrect) mean(AUC_later_incorrect)],[std(AUC_1_correct) std(AUC_later_correct) std(AUC_1_incorrect) std(AUC_later_incorrect)])
hold on
box off
ylim([0 1])
xlim([0.8 4.2])
plot([0.8 4.2],[0.5 0.5])

%% ANOVA test for figure 2D

%% unbalanced
y=[AUC_1_correct,AUC_1_incorrect,AUC_later_correct,AUC_later_incorrect]';
correct=[ones(size(AUC_1_correct)),zeros(size(AUC_1_incorrect)),ones(size(AUC_later_correct)),zeros(size(AUC_later_incorrect))]';
first_later=[ones(size(AUC_1_correct)),ones(size(AUC_1_incorrect)),ones(size(AUC_later_correct))*2,ones(size(AUC_later_incorrect))*2]';

p = anovan(y,{correct,first_later});
p(1)
p(2)
pause
% plot(AUC_later_correct,AUC_1_correct,'ob')
% plot(mean(AUC_later_correct),mean(AUC_1_correct),'.b')
% errorbar(mean(AUC_later_correct),mean(AUC_1_correct),std(AUC_later_correct),'horizontal','c')
% errorbar(mean(AUC_later_correct),mean(AUC_1_correct),std(AUC_1_correct),'c')
% 
% plot(AUC_later_incorrect,AUC_1_incorrect,'or')
% plot(mean(AUC_later_incorrect),mean(AUC_1_incorrect),'.r')
% errorbar(mean(AUC_later_incorrect),mean(AUC_1_incorrect),std(AUC_later_incorrect),'horizontal','m')
% errorbar(mean(AUC_later_incorrect),mean(AUC_1_incorrect),std(AUC_1_incorrect),'m')
% box off
% xlim([0 1])
% ylim([0 1])
% plot([0 1],[0 1],'k')
% plot([0 0.5],[0.5 0.5],'--k')
% plot([0.5 0.5],[0 0.5],'--k')

xlabel('AUC later touches')
ylabel('AUC first touch')
hold on
% 
% plot(AUC_later_incorrect(3),AUC_1_incorrect(3),'sm')
% plot(AUC_later_correct(3),AUC_1_correct(3),'sc')

%subplot(2,2,4)
%bar([AUC_1;AUC_later2]')
%box off
%xlabel('Session number')
%ylabel('AUC')
%hold on

cd ..
end

function present=onthelist(txt,unit)

for i=1:size(txt,1)
    Index=strcmp([txt{i,:}],unit);
    if Index==0
        present=0;
    else
        present=1;
        return
    end
end
end