%% The aim of this script is to reproduce the figures of the paper 
%% This script should be run in the folder "Data", which contains the subfolders (one for each session) with the behavioural and neural data.
close all
clear all

%Figure 1
variability_behaviour

% Figure 2

% Figure 3

% Figure 4
[fig4,change_explained_by_adap]=all_Tcurves;
gradual_adaptation(fig4)

% Figure 5
