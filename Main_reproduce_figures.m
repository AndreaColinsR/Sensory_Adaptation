%% The aim of this script is to reproduce the figures of the paper 
%% This script should be run in the folder "Data", which contains the subfolders (one for each session) with the behavioural and neural data.
close all
clear all

%Figure 1
variability_behaviour

% Figure 2
[FR_touch,FR_touch_correct,FR_touch_incorrect]=all_Tcurves_correct_incorrect_population;
fig2=figure;


% Figure 3

% Figure 4
[fig4,change_explained_by_adap,FR_unit]=all_Tcurves;
gradual_adaptation(fig4)


figure(fig2)

subplot(4,3,3)
plot(1:4,FR_unit','Color',[0.5 0.5 0.5])
hold on 
% errorbar(1:4,mean(FR_touch,1),std(FR_touch,[],1),'k')
errorbar(1:4,mean(FR_touch_correct,1),std(FR_touch_correct,[],1),'g')
errorbar(1:4,mean(FR_touch_incorrect,1),std(FR_touch_incorrect,[],1),'r')
xlabel('Touch number')
xticks([1 2 3 4])
ylabel('FR [spikes/touch]')
xlim([0.5 4.5])
box off

% Figure 5
