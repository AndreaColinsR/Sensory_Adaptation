function recovery_test(fig2)
folder=[{'Glu32_19092017H'};{'Glu32_21092017H'};{'Glu43_22122017H'};{'Glu35_10112017H'};{'Glu35_13112017H_1'};{'Glu35_13112017H_2'}];

[~,txt]=xlsread('Adaptation_units_list.xlsx',1);

nunits=size(txt,1);

AUC_1_correct=nan(1,nunits);
AUC_later_correct=nan(1,nunits);
AUC_12_correct=nan(1,nunits);
AUC_later2_correct=nan(1,nunits);


AUC_1_incorrect=nan(1,nunits);
AUC_later_incorrect=nan(1,nunits);
AUC_12_incorrect=nan(1,nunits);
AUC_later2_incorrect=nan(1,nunits);

AUC_1=nan(1,nunits);
AUC_later=nan(1,nunits);
AUC_12=nan(1,nunits);
AUC_later2=nan(1,nunits);
counter=1;
figure(fig2)
for f=1:size(folder,1)
    
    
    
    cd(folder{f})
    info=xlsread('ledtrials.xlsx');
    go_trials=find(info(:,7)==2);
    idx_correct=info(go_trials,6)==1;
    idx_incorrect=info(go_trials,6)==2;
    
    ff=dir('*neural_data_S*');
    load(ff(1).name)
    
    for unit=1:size(ff,1)
        
        if onthelist(txt,[folder{f} ff(unit).name])
            load(ff(unit).name)
            total_psth=psth;
            
            if counter==3
                do_extra_plot=1;
            else
                do_extra_plot=0;
            end
            
            
            [hit_correct,hit3_correct,fa_correct,fa2_correct,fa3_correct]=Detection_test(total_psth(idx_correct,:),idx_correct);
            
            AUC_1_correct(counter)=trapz(fliplr([ 1 fa_correct]),fliplr([ 1 hit_correct]));
            AUC_later_correct(counter)=trapz(fliplr([ 1 fa3_correct]),fliplr([ 1 hit3_correct]));
            AUC_12_correct(counter)=trapz(fliplr([ 1 fa2_correct]),fliplr([ 1 hit_correct]));
            AUC_later2_correct(counter)=trapz(fliplr([ 1 fa3_correct]),fliplr([ 1 hit3_correct]));
            
            if do_extra_plot
                subplot(4,3,9)
                plot([ 1 fa_correct] , [1 hit_correct],'-r')
                hold on
                plot([ 1 fa3_correct] , [1 hit3_correct],'-k')
                plot([0 1],[0 1],'b')
                box off
                
            end
            
            
            [hit_incorrect,hit3_incorrect,fa_incorrect,fa2_incorrect,fa3_incorrect]=Detection_test(total_psth(idx_incorrect,:),idx_incorrect);
            
            AUC_1_incorrect(counter)=trapz(fliplr([ 1 fa_incorrect]),fliplr([ 1 hit_incorrect]));
            AUC_later_incorrect(counter)=trapz(fliplr([ 1 fa3_incorrect]),fliplr([ 1 hit3_incorrect]));
            AUC_12_incorrect(counter)=trapz(fliplr([ 1 fa2_incorrect]),fliplr([ 1 hit_incorrect]));
            AUC_later2_incorrect(counter)=trapz(fliplr([ 1 fa3_incorrect]),fliplr([ 1 hit3_incorrect]));
            
            if do_extra_plot
                plot([ 1 fa_incorrect] , [1 hit_incorrect],'--r')
                plot([ 1 fa3_incorrect] , [1 hit3_incorrect],'--k')
                hold off
            end
            
            [hit,hit3,fa,fa2,fa3]=Detection_test(total_psth);
            AUC_1(counter)=trapz(fliplr([ 1 fa]),fliplr([ 1 hit]));
            AUC_later(counter)=trapz(fliplr([ 1 fa3]),fliplr([ 1 hit3]));
            AUC_12(counter)=trapz(fliplr([ 1 fa2]),fliplr([ 1 hit]));
            AUC_later2(counter)=trapz(fliplr([ 1 fa3]),fliplr([ 1 hit3]));
            
            counter=counter+1;
            
            clear total_psth
        end
    end
    cd ..
end
subplot(4,3,9)
xlabel('False Alarm')
ylabel('Hit')

subplot(4,3,12)
errorbar([0.8 1.8],[mean(AUC_1_correct) mean(AUC_1_incorrect)],[std(AUC_1_correct) std(AUC_1_incorrect)],'.r')
hold on
errorbar([1.2 2.2],[ mean(AUC_later_correct) mean(AUC_later_incorrect)],[std(AUC_later_correct) std(AUC_later_incorrect)],'.k')
box off
ylim([0 1])
xlim([0.6 2.4])
plot([0 2.5],[0.5 0.5],'Color',[0.5 0.5 0.5])
xticks([1 2])
xticklabels({'Hit','Miss'})
ylabel('AUC')
hold on

%% ANOVA test for figure 2D

%% unbalanced
y=[AUC_1_correct,AUC_1_incorrect,AUC_later_correct,AUC_later_incorrect]';
correct=[ones(size(AUC_1_correct)),zeros(size(AUC_1_incorrect)),ones(size(AUC_later_correct)),zeros(size(AUC_later_incorrect))]';
first_later=[ones(size(AUC_1_correct)),ones(size(AUC_1_incorrect)),ones(size(AUC_later_correct))*2,ones(size(AUC_later_incorrect))*2]';

p = anovan(y,{correct,first_later});
p(1)
p(2)




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