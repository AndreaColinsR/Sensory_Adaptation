% %% The aim of this script is to reproduce the figures of the paper 
% %% This script should be run in the folder "Data", which contains one matlab file per recording with the behavioural and neural data.
% 
close all
clearvars
rng('default') % for reproducibility
% tic
%% Figure 1
fig1=figure;
variability_behaviour(fig1)

%% Figure 2
fig2=figure;
[FR_touch,FR_touch_correct,FR_touch_incorrect]=all_Tcurves_correct_incorrect_population;
recovery_test(fig2)


%% Figure 3
fig3=figure;
strength_of_whisker_touches(fig3)
distance_base_all_sessions(fig3)

%% Figure 4

fig4=figure;
[change_explained_by_adap,FR_unit]=all_Tcurves(fig4);
gradual_adaptation(fig4)


% % Figure 2 panel C
figure(fig2)
subplot(4,3,3)
plot(1:4,FR_unit','Color',[0.5 0.5 0.5])
hold on 
errorbar(1:4,mean(FR_touch,1),std(FR_touch,[],1),'k')
errorbar(1:4,mean(FR_touch_correct,1),std(FR_touch_correct,[],1),'g')
errorbar(1:4,mean(FR_touch_incorrect,1),std(FR_touch_incorrect,[],1),'r')
xlabel('Touch number')
xticks([1 2 3 4])
ylabel('FR [spikes/touch]')
xlim([0.5 4.5])
box off

%% Figure 5
tic
fig5=figure;
cross_whisker_adaptation(fig5)
toc