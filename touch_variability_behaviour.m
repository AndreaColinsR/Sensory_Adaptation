function [nt_per_trial,inter_time_touch,length_touch]=touch_variability_behaviour(Data,selected_trials,do_plot)
%% touch_variability_behaviour calculates the most important touch parameters for the selected trials in Data
% INPUTS
% Data: data structure containing the behavioural and neural data
%
% selected_trials: binary array indicating the trials selected for this
% analysis
%
% do_plot: 1-plot example trial (dk and touch epochs)
%          0-do not plot
%
% OUTPUTS
%
% nt_per_trial: number of touches per trial
%
% inter_time_touch: time between the start of two consecutive touches [ms]
%
% length of touch: time that the whisker is in contact with the pole [ms]
% 

touches_whisker=Data.touch_per_whisker(selected_trials,:,:);

ntrials=size(touches_whisker,1);
k_c1=Data.deltak_w1(selected_trials,:);
k_c2=Data.deltak_w2(selected_trials,:);

a_c1=Data.azimuth_w1(selected_trials,:);
a_c2=Data.azimuth_w2(selected_trials,:);

touches_m=round(ones(ntrials,size(k_c1,2),3))*255;
counter=1;
nt_per_trial=zeros(ntrials,1);
nt_per_trial_c1=zeros(ntrials,1);
nt_per_trial_c2=zeros(ntrials,1);
counter2=1;
for i=1:ntrials
    
    
    %% for c1
    touches_c1=touches_whisker(i,:,1);
    not_touch=~touches_c1;
    touches_m(i,:,1)=not_touch.*255;
    touches_m(i,:,2)=not_touch.*255;
    clear not_touch touches
    
    touch_onset=find(diff(touches_c1)>0);
    if ~isempty(touch_onset)
        nt_per_trial(i)=numel(touch_onset);
        nt_per_trial_c1(i)=nt_per_trial(i);
    end
    for j=1:numel(touch_onset)
        
        tmp=find(diff(touches_c1(touch_onset(j):end))<0,1,'First');
        if isempty(tmp)
            tmp=3488-touch_onset(j);
            length_touch(counter)=tmp;
            
        else
            length_touch(counter)=tmp;
            
        end
        counter=counter+1;
    end
    
    %% for c2
    
    touches_c2=touches_whisker(i,:,2);
    not_touch=~touches_c2;
    touches_m(i,:,1)=(not_touch & touches_m(i,:,1)).*255;
    touches_m(i,:,3)=not_touch.*255;
    
    touch_onset=find(diff(touches_c2)>0);
    if ~isempty(touch_onset)
        nt_per_trial(i)= nt_per_trial(i)+numel(touch_onset);
        nt_per_trial_c2(i)=numel(touch_onset);
    end
    
    for j=1:numel(touch_onset)
        
        tmp=find(diff(touches_c2(touch_onset(j):end))<0,1,'First');
        if isempty(tmp)
            tmp=3488-touch_onset(j);
            length_touch(counter)=tmp;
            
        else
            length_touch(counter)=tmp;
            
        end
        counter=counter+1;
    end
    
    clear not_touch touches
    
    %% inter-touch interval
    touches_all=sum(touches_whisker(i,:,:),3);
    
    touch_offset=find(diff(touches_all)<0);
    
    for j=1:numel(touch_offset)
        if touch_offset(j)<3488 && (touches_all(touch_offset(j)+1)==0)
            new_onset=find(diff(touches_all(touch_offset(j):end))>0,1,'First');
            
            if ~isempty(new_onset)
                inter_time_touch(counter2)=new_onset;
                counter2=counter2+1;
                
            end
            
        end
        
    end
    
    if do_plot && i==9
        
    subplot(4,2,1)
    image([0 size(touches_m,2)],[-0.03 0.02],touches_m(i,:,:),'AlphaData',0.2)
    hold on 
    plot(k_c1(i,:),'b')
    plot(k_c2(i,:),'g')
    ylim([-0.03 0.02])
    set(gca,'YDir','normal')
    legend('C1','C2')
    hold off
    box off

    subplot(4,2,3)
    hold on
    image([0 size(touches_m,2)],[0 100],touches_m(i,:,:),'AlphaData',0.2)
    plot(a_c1(i,:),'b')
    plot(a_c2(i,:),'g')
    set(gca,'YDir','normal')
    legend('C1','C2')
    ylim([20 100])
    xlim([0 3488])
    hold off
    box off

    end
end

if do_plot

    subplot(4,2,[5 7])
    image(touches_m,'AlphaData',0.2)
    xlim([0 3488])
    xlabel('Time from trial start [ms]')
    ylabel('Trials')
    set(gca, 'YDir','normal')
    ylim([0 84])
    box off
end

% figure
% subplot(4,1,2)
% histogram(length_touch(length_touch<=350),[0:10:350],'Normalization','Probability')
% n_touches=numel(length_touch);
% over_400=100*sum(length_touch>350)/n_touches
% hold on
% xlabel('Touch length')
% ylabel('Fraction of touches')
% hold off
% box off
%
% subplot(4,1,1)
% histogram(nt_per_trial(nt_per_trial>0),'Normalization','Probability')
% hold on
% xlabel('Touches per trial')
% ylabel('Fraction of trials')
% box off
%
% subplot(4,1,4)
% bar([sum(nt_per_trial_c1) sum(nt_per_trial_c2)]/sum(nt_per_trial_c1+nt_per_trial_c2))
% hold on
% xlabel('Whisker')
% ylabel('Fraction of touches')
% box off
%
% subplot(4,1,3)
% histogram(inter_time_touch,[0:10:350],'Normalization','Probability')
% hold on
% over_400=100*sum(inter_time_touch>350)/numel(inter_time_touch)
% xlabel('Inter-touch interval')
% ylabel('Fraction of touches')
% box off



% hold off
%
% subplot(3,1,3)
% h1=histogram(abs(slip_detach),[0:1:200]);
% 100*sum(h1.Values(1:5))/(nr_files-t_no_slips)
% hold on
% xlabel('Time between the first slip and the first detach [ms]')
% ylabel('Number of trials')
% hold off
%make an histogram of the diference in time first slip-first touch
end