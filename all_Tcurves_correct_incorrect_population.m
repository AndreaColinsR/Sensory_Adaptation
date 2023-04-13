function [FR_touch,FR_touch_correct,FR_touch_incorrect]=all_Tcurves_correct_incorrect_population
counter=0;
av_unit=0;
ff=dir('*.mat*');
Nsessions=size(ff,1);

percentile_1=0.95;
percentile_2=0.95;

FR_correct_1_all=[];
FR_incorrect_1_all=[];
FR_correct_later_all=[];
FR_incorrect_later_all=[];

Nsamples_tuning_curve=nan(Nsessions,2);

%% Define population FR per mouse
mouse_id=zeros(Nsessions,1);

for f=1:Nsessions
mouse_id(f)=str2double(ff(f).name(4:5));
end

mouse_number=unique(mouse_id);
FR_mouse=cell(numel(mouse_number),1);
adapt_mouse=cell(numel(mouse_number),1);


for f=1:Nsessions
    
    load(ff(f).name,'Data')
    this_mouse=str2double(ff(f).name(4:5))==mouse_number;

    %% select correct and incorrect trials
    touches_matrix_correct=Data.touch(Data.correct_trials,:);
    
    touches_whisker_correct=Data.touch_per_whisker(Data.correct_trials,:,:);
    
    k_c1_correct=Data.deltak_w1(Data.correct_trials,:);
    k_c2_correct=Data.deltak_w2(Data.correct_trials,:);
    touches_matrix_incorrect=Data.touch(Data.incorrect_trials,:);
    touches_whisker_incorrect=Data.touch_per_whisker(Data.incorrect_trials,:,:);
    k_c1_incorrect=Data.deltak_w1(Data.incorrect_trials,:);
    k_c2_incorrect=Data.deltak_w2(Data.incorrect_trials,:);
    total_psth=zeros(size(Data.deltak_w1));
    total_psth_correct=zeros(size(k_c1_correct));
    total_psth_incorrect=zeros(size(k_c1_incorrect));
    
    for i_unit=1:size(Data.unit,2)
        
            total_psth=Data.unit(i_unit).spikes+total_psth;
            total_psth_correct=Data.unit(i_unit).spikes(Data.correct_trials,:)+total_psth_correct;
            total_psth_incorrect=Data.unit(i_unit).spikes(Data.incorrect_trials,:)+total_psth_incorrect;
            
            
            delete=size(total_psth,2)-size(Data.deltak_w2,2);
            Data.deltak_w1=[zeros(size(Data.deltak_w1,1),delete),Data.deltak_w1];
            Data.deltak_w2=[zeros(size(Data.deltak_w2,1),delete),Data.deltak_w2];
            clear deltak touch_idx FR

    end
    
    [x1,p1,x2,p2,x1_norm,p1_norm,x2_norm,p2_norm,x_all,p_all,deltak,touch_idx,FR,~,FR_prev,whisker]=Tcurve_touches_dk(Data.touch,Data.touch_per_whisker,total_psth,Data.deltak_w1,Data.deltak_w2,0,4,0,percentile_1,percentile_2);
    [x1_correct,p1_correct,x2_correct,p2_correct,~,~,~,~,~,~,~,touch_idx_correct,FR_correct,~,~,~,spikes_correct]=Tcurve_touches_dk(touches_matrix_correct,touches_whisker_correct,total_psth_correct,k_c1_correct,k_c2_correct,0,4,0,percentile_1,percentile_2);
    [x1_incorrect,p1_incorrect,x2_incorrect,p2_incorrect,~,~,~,~,~,~,~,touch_idx_incorrect,FR_incorrect,~,~,~,spikes_incorrect]=Tcurve_touches_dk(touches_matrix_incorrect,touches_whisker_incorrect,total_psth_incorrect,k_c1_incorrect,k_c2_incorrect,0,4,0,percentile_1,percentile_2);
    
    deltak=abs(deltak);
    
    if strcmp('Glu32_19092017H.mat',ff(f).name)
    raster_plot(spikes_correct,touch_idx_correct,[1 2 4 5])
    raster_plot(spikes_incorrect,touch_idx_incorrect,[7 8 10 11])
    end
    
    
    
    
    [lambda1,~]=Tcurve_touches_dk_cross_val(Data.touch,Data.touch_per_whisker,total_psth,Data.deltak_w1,Data.deltak_w2,0,4);
    
    if ~isnan(p1_incorrect)
        av_unit=av_unit+1;
        
%         spikes_incorrect_1(av_unit)=sum(FR_correct(touch_idx_correct==1));
%         coeff1(av_unit,:)=p1;
%         coeff2(av_unit,:)=p2;
        
        
%         f1_correct = polyval(p1_correct,x1_correct);
%         f2_correct = polyval(p2_correct,x1_correct);
%         f1_incorrect = polyval(p1_incorrect,x1_incorrect);
%         f2_incorrect = polyval(p2_incorrect,x1_incorrect);
%         pos_slope(av_unit)=p_all(1);
%         
%         
%         slope_first(av_unit,:)=[p1_correct(1) p1_incorrect(1)];
%         slope_later(av_unit,:)=[p2_correct(1) p2_incorrect(1)];
%         intercept_first(av_unit,:)=[p1_correct(2), p1_incorrect(2)];
%         intercept_later(av_unit,:)=[p2_correct(2), p2_incorrect(2)];
        
        Nsamples_tuning_curve(f,:)=[sum(touch_idx==1) sum(touch_idx>1)];
        
        FR1_all=simspikes(deltak(touch_idx==1)',p2,1); %first touch with mixed tunning cuve
        FR2=simspikes(deltak(touch_idx>1)',p2,1); %later touches with corresponding tuning curve
        FR2_all=simspikes(deltak(touch_idx>1)',p_all,1); %later touches with mixed tuning curve
        
        
        %ratioFR(av_unit,:)=[mean(FR(touch_idx==1)'./lambda1) 1/nanmean(FR(touch_idx==1)'./FR1_all(:)')];
        % response to first touch is X% higher than predicted
        %ratioFRb(av_unit,:)=nanmean(FR(touch_idx==1)'./FR1_all(:)');
        
        %% test if FR from first touch is higher than prediction
        [hFR(av_unit),pvalFR(av_unit)] = ttest(FR(touch_idx==1)',FR1_all(:)','Tail','right');
        
        
        %error1(av_unit)=mean(abs(FR(touch_idx==1)-lambda1(:)));
        %error2(av_unit)=mean(abs(FR(touch_idx==1)-FR1_all(:)));
        
        %ratioFR2(av_unit,:)=[1 mean(FR2(:))/mean(FR2_all(:))];
        
        
        FR_correct_incorrect_1(av_unit,:)=[mean(FR_correct(touch_idx_correct==1)) mean(FR_incorrect(touch_idx_incorrect==1))];
        FR_correct_incorrect_later(av_unit,:)=[mean(FR_correct(touch_idx_correct>1)) mean(FR_incorrect(touch_idx_incorrect>1))];
        
        %all trials
        FR_correct_1_all=[FR_correct_1_all;FR_correct(touch_idx_correct==1)];
        FR_incorrect_1_all=[FR_incorrect_1_all;FR_incorrect(touch_idx_incorrect==1)];
        FR_correct_later_all=[FR_correct_later_all;FR_correct(touch_idx_correct>1)];
        FR_incorrect_later_all=[FR_incorrect_later_all;FR_incorrect(touch_idx_incorrect>1)];
        
        
        hw(av_unit)=ttest2(FR((whisker==1)),FR((whisker==2)));
        
        [h_atte(av_unit),pval_atte(av_unit)]=ttest2(FR(touch_idx==1),FR(touch_idx>1));
        
            change_explained_by_adap(av_unit)=min([(median(FR(touch_idx==1)-FR1_all(:)))/(median(FR(touch_idx==1))-median(FR(touch_idx>1))) 1e6]);
            
            adapt_mouse{this_mouse}= [adapt_mouse{this_mouse};change_explained_by_adap(av_unit)];
        
        FR_touch(av_unit,:)=[mean(FR(touch_idx==1)) mean(FR(touch_idx==2)) mean(FR(touch_idx==3)) mean(FR(touch_idx>3))]/mean(FR(touch_idx==1));
        FR_touch_correct(av_unit,:)=[mean(FR_correct(touch_idx_correct==1)) mean(FR_correct(touch_idx_correct==2)) mean(FR_correct(touch_idx_correct==3)) mean(FR_correct(touch_idx_correct>3))]/mean(FR_correct(touch_idx_correct==1));
        FR_touch_incorrect(av_unit,:)=[mean(FR_incorrect(touch_idx_incorrect==1)) mean(FR_incorrect(touch_idx_incorrect==2)) mean(FR_incorrect(touch_idx_incorrect==3)) mean(FR_incorrect(touch_idx_incorrect>3))]/mean(FR_incorrect(touch_idx_incorrect==1));
        
        FR_mouse{this_mouse}= [FR_mouse{this_mouse};FR_touch(av_unit,:)];
        
        if p2(1)<p1(1)
            counter=counter+1;
        end
        
        clear FR FR_prev FR1_all FR2_all FR2
    end
    
    clear touch_idx
    
end


subplot(4,3,6)
hold on
% hit
errorbar(0.8,mean(FR_correct_1_all),std(FR_correct_1_all) ,'-.m')
errorbar(1.2,mean(FR_correct_later_all),std(FR_correct_later_all),'-.b')

% miss
errorbar(1.8,mean(FR_incorrect_1_all),std(FR_incorrect_1_all),'-.m')
errorbar(2.2,mean(FR_incorrect_later_all), std(FR_incorrect_later_all),'-.b')

ylabel('Firing rate [spikes/touch]')
box off
xlim([0.6 2.4])
xticks([1 2])
xticklabels({'Hit','Miss'})
plot([0.8 1.8],FR_correct_incorrect_1,'.r')
plot([1.2 2.2],FR_correct_incorrect_later,'.k')


% [h_FR_1,pFR_1] = ttest2(FR_correct_1_all,(FR_incorrect_1_all));
% [h_FR_later,pFR_later] = ttest2(FR_correct_later_all,FR_incorrect_later_all);

%
% [h_FR_correct,pFR_correct] = ttest2(FR_correct_1_all,FR_correct_later_all);
% [h_FR_incorrect,pFR_incorrect] = ttest2(FR_incorrect_1_all,FR_incorrect_later_all);
%pause

%% ANOVA Fig 2C (FRs)
%% unbalanced
correct=[ones(size(FR_correct_1_all));zeros(size(FR_incorrect_1_all));ones(size(FR_correct_later_all));zeros(size(FR_incorrect_later_all))];
first_later=[ones(size(FR_correct_1_all));ones(size(FR_incorrect_1_all));ones(size(FR_correct_later_all))*2;ones(size(FR_incorrect_later_all))*2];
y=[FR_correct_1_all;FR_incorrect_1_all;FR_correct_later_all;FR_incorrect_later_all];


p = anovan(y,{correct,first_later},'display','off');

disp('-----------------------------------')
disp('Unbalanced anova First vs later responses (FR)')
disp(['p-value Effect hit vs miss = ' num2str(p(1))])
disp(['p-value Effect first vs later = ' num2str(p(2))])

%% balanced
nsample=size(FR_incorrect_1_all,1);
balanced_anova(FR_correct_1_all,FR_incorrect_1_all,FR_correct_later_all,FR_incorrect_later_all,nsample)


%%
% fraction between Responses in miss trials and hit trials.
disp(' ')
disp('-----------------------------------')
miss_hit_ratio=mean([FR_incorrect_1_all;FR_incorrect_later_all])/mean([FR_correct_1_all;FR_correct_later_all]);
disp(['Ration responses miss/hit = ' num2str(miss_hit_ratio)])

disp('-----------------------------------')

%%%%% ttest FR
%errorbar([1 2],mean(FR_correct_incorrect_1),std(FR_correct_incorrect_1),'-m')
%errorbar([1 2],mean(FR_correct_incorrect_later),std(FR_correct_incorrect_later),'-b')
%[h_FR_1,pFR_1] = ttest2(FR_correct_incorrect_1(:,1),FR_correct_incorrect_1(:,2));
%[h_FR_later,pFR_later] = ttest2(FR_correct_incorrect_later(:,1),FR_correct_incorrect_later(:,2));



% [h1,p1] = ttest(FR_touch(:,2),1,'Tail','left');
% [h2,p2] = ttest2(FR_touch(:,2),FR_touch(:,3),'Tail','right');
% [h3,p3] = ttest2(FR_touch(:,3),FR_touch(:,4),'Tail','right');


disp('')
disp('-----------------------------------------')
disp('Change due to adaptation (alpha) whole population')
disp(num2str(mean(change_explained_by_adap)))


disp('')
disp('')
disp('-----------------------------------------')
disp('Change due to adaptation (alpha) per mouse')

for im=1:numel(mouse_number)
    
disp(['Mouse ' num2str(im) ' =' num2str(mean(adapt_mouse{im}))])

end


%% Report decrease in FR per mouse
disp('')
disp('')
disp('-----------------------------------------')
disp('Normalised FR of second touch (FR2/FR1)')
for im=1:numel(mouse_number)
    
disp(['Mouse ' num2str(im) ' =' num2str(mean(FR_mouse{im}(:,2),1))])

end


%% Number of touches for each tuning curve
disp('')
disp('')
disp('-----------------------------------------')
disp('Average number of touches for each tuning curve')
disp([' First touch = ' num2str(mean(Nsamples_tuning_curve(:,1))) '  Later touches = ' num2str(mean(Nsamples_tuning_curve(:,2)))])

end

function FR=simspikes(x,p,Nsim)
lambda=p(1)*x+p(2);
lambda=repmat(lambda,1,Nsim);
FR = lambda;
end
