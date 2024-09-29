function recovery_test(fig2)
% this function test if later touches are less discriminable from the noise
% (Fig 2)
ff=dir('*.mat*');
nunits=33;

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
for f=1:size(ff,1)
    
    load(ff(f).name,'Data')
    
    
    for i_unit=1:size(Data.unit,2)
        
        %if onthelist(txt,[folder{f} ff(unit).name])
        %    load(ff(i_unit).name)
            total_psth=Data.unit(i_unit).spikes;
            
            if counter==3
                do_extra_plot=1;
            else
                do_extra_plot=0;
            end
            
            
            [hit_correct,hit3_correct,fa_correct,fa2_correct,fa3_correct]=Detection_test(Data,total_psth(Data.correct_trials,:),Data.correct_trials);
            
            AUC_1_correct(counter)=trapz(fliplr([ 1 fa_correct]),fliplr([ 1 hit_correct]));
            AUC_later_correct(counter)=trapz(fliplr([ 1 fa3_correct]),fliplr([ 1 hit3_correct]));
            AUC_12_correct(counter)=trapz(fliplr([ 1 fa2_correct]),fliplr([ 1 hit_correct]));
            AUC_later2_correct(counter)=trapz(fliplr([ 1 fa3_correct]),fliplr([ 1 hit3_correct]));
            
            if do_extra_plot
                subplot(4,6,18)
                plot([ 1 fa_correct] , [1 hit_correct],'-r')
                hold on
                plot([ 1 fa3_correct] , [1 hit3_correct],'-k')
                plot([0 1],[0 1],'b')
                box off
                
            end
            
            
            [hit_incorrect,hit3_incorrect,fa_incorrect,fa2_incorrect,fa3_incorrect]=Detection_test(Data,total_psth(Data.incorrect_trials,:),Data.incorrect_trials);
            
            AUC_1_incorrect(counter)=trapz(fliplr([ 1 fa_incorrect]),fliplr([ 1 hit_incorrect]));
            AUC_later_incorrect(counter)=trapz(fliplr([ 1 fa3_incorrect]),fliplr([ 1 hit3_incorrect]));
            AUC_12_incorrect(counter)=trapz(fliplr([ 1 fa2_incorrect]),fliplr([ 1 hit_incorrect]));
            AUC_later2_incorrect(counter)=trapz(fliplr([ 1 fa3_incorrect]),fliplr([ 1 hit3_incorrect]));
            
            if do_extra_plot
                plot([ 1 fa_incorrect] , [1 hit_incorrect],'--r')
                plot([ 1 fa3_incorrect] , [1 hit3_incorrect],'--k')
                hold off
            end
            
            [hit,hit3,fa,fa2,fa3]=Detection_test(Data,total_psth);
            AUC_1(counter)=trapz(fliplr([ 1 fa]),fliplr([ 1 hit]));
            AUC_later(counter)=trapz(fliplr([ 1 fa3]),fliplr([ 1 hit3]));
            AUC_12(counter)=trapz(fliplr([ 1 fa2]),fliplr([ 1 hit]));
            AUC_later2(counter)=trapz(fliplr([ 1 fa3]),fliplr([ 1 hit3]));
            
            counter=counter+1;
            
            clear total_psth
        %end
    end
    %cd ..
end
subplot(4,6,18)
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

y=[AUC_1_correct,AUC_1_incorrect,AUC_later_correct,AUC_later_incorrect]';
correct=[ones(size(AUC_1_correct)),zeros(size(AUC_1_incorrect)),ones(size(AUC_later_correct)),zeros(size(AUC_later_incorrect))]';
first_later=[ones(size(AUC_1_correct)),ones(size(AUC_1_incorrect)),ones(size(AUC_later_correct))*2,ones(size(AUC_later_incorrect))*2]';

p = anovan(y,{correct,first_later},'display','off');

disp('-----------------------------------')
disp('Anova First vs later touch detection (AUC)')
disp(['p-value Effect hit vs miss = ' num2str(p(1))])
disp(['p-value Effect first vs later = ' num2str(p(2))])


end