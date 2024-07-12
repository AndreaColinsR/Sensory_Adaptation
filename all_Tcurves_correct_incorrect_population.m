function [FR_touch,FR_touch_correct,FR_touch_incorrect]=all_Tcurves_correct_incorrect_population
counter=0;
ff=dir('*.mat*');
Nsessions=size(ff,1);

percentile=0.95;
NpointsTcurve=4;

Nsamples_tuning_curve=nan(Nsessions,2);
FR_correct_incorrect_1=nan(Nsessions,2);
FR_correct_incorrect_later=nan(Nsessions,2);
change_explained_by_adap=nan(Nsessions,1);  
pval_atte=nan(Nsessions,1);  
pvalFR=nan(Nsessions,1);
FR_touch=nan(Nsessions,4);
FR_touch_correct=nan(Nsessions,4);
FR_touch_incorrect=nan(Nsessions,4);

FR_all.correct_1=[];
FR_all.incorrect_1=[];
FR_all.correct_later=[];
FR_all.incorrect_later=[];

        

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
    
    [~,p1,~,p2,~,~,deltak,touch_idx,FR,~,~,whisker]=Tcurve_touches_dk(Data.touch,Data.touch_per_whisker,total_psth,Data.deltak_w1,Data.deltak_w2,0,NpointsTcurve,0,percentile,percentile);
    [~,~,~,~,~,~,~,touch_idx_correct,FR_correct,~,~,~,spikes_correct]=Tcurve_touches_dk(touches_matrix_correct,touches_whisker_correct,total_psth_correct,k_c1_correct,k_c2_correct,0,NpointsTcurve,0,percentile,percentile);
    [~,p1_incorrect,~,~,~,~,~,touch_idx_incorrect,FR_incorrect,~,~,~,spikes_incorrect]=Tcurve_touches_dk(touches_matrix_incorrect,touches_whisker_incorrect,total_psth_incorrect,k_c1_incorrect,k_c2_incorrect,0,NpointsTcurve,0,percentile,percentile);
    deltak=abs(deltak);


    %% adding latencies analyses
    [Lat,touch_lat,prev_FR]=Latencies(Data.touch,Data.touch_per_whisker,total_psth,1,2);

    %% Are initial latencies different for different touch numbers?
    [~,pvalLat]=ttest2(Lat(touch_lat==1),Lat(touch_lat>=3));
    disp(['p-value initial latencies (first touch vs late touch) = ' num2str(pvalLat,'%.4f')])

    [~,pvalFRLat]=ttest2(prev_FR(touch_lat==1),prev_FR(touch_lat>=3));

    %disp([' corr latency vs FR pre touch =' num2str(corr(Lat',prev_FR'))])

    %errorbar([1 4],[mean(prev_FR(touch_lat==1)),mean(prev_FR(touch_lat>=4))],[std(prev_FR(touch_lat==1)),std(prev_FR(touch_lat>=4))])
    %disp(['p-value FR prev (first touch vs late touch) = ' num2str(pvalFRLat,'%.4f')])
    
    %% Latencies early vs late trials
    %Ntrials=size(Data.touch,1);
    %Early_trials=1:floor(Ntrials/2);
    %Late_trials=floor(Ntrials/2)+1:Ntrials;
    %[Lat_early,touch_lat_early]=Latencies(Data.touch(Early_trials,:),Data.touch_per_whisker(Early_trials,:,:),total_psth(Early_trials,:),0,2);
    
    %[Lat_late,touch_lat_late]=Latencies(Data.touch(Late_trials,:),Data.touch_per_whisker(Late_trials,:,:),total_psth(Late_trials,:),0,2);
    %[~,pvalLat]=ttest2(Lat_early,Lat_late);
    %disp(['p-value early vs late trials latencies = ' num2str(pvalLat,'%.4f')])
    
    %% plot example raster plot
    if strcmp('Glu32_19092017H.mat',ff(f).name)
        raster_plot(spikes_correct,touch_idx_correct,[1 2 4 5])
        raster_plot(spikes_incorrect,touch_idx_incorrect,[7 8 10 11])
    end
    
    
    if ~isnan(p1_incorrect)
        
        %% testing whisker preference
        % does this neuron show a preference for a whisker
%         hw(av_unit)=ttest2(FR((whisker==1 & touch_idx>4)),FR((whisker==2 & touch_idx>4)));
%         
%         if hw(av_unit)
%             disp(['Preference session =' num2str(av_unit)])
%             % choose the prefered whisker
%             %[~,pref_w]=max([mean(FR((whisker==1 & touch_idx>4))),mean(FR((whisker==2 & touch_idx>4)))])
%             
%             % choose the whisker that touches more times
%             [~,pref_w]=max([sum(whisker==1 & touch_idx>4),sum((whisker==2 & touch_idx>4))]);
%                 
%             
%             selected_touch=touch_idx==1 & whisker==pref_w;
%             nsamples_1=sum(selected_touch)
%             [~,~,p1]=plotTCurve(deltak(selected_touch),FR(selected_touch),'',NpointsTcurve,percentile);
%             
%             selected_touch=touch_idx>1 & whisker==pref_w;
%             nsamples_2=sum(selected_touch)
%             [~,~,p2]=plotTCurve(deltak(selected_touch),FR(selected_touch),'',NpointsTcurve,percentile);
%             
%             FR=FR(whisker==pref_w);
%             deltak=deltak(whisker==pref_w);
%             touch_idx=touch_idx(whisker==pref_w);
%         end
        
        
        
        Nsamples_tuning_curve(f,:)=[sum(touch_idx==1) sum(touch_idx>1)];
        
        FR1_all=simspikes(deltak(touch_idx==1)',p2,1); %first touch with mixed tunning cuve
        %FR2=simspikes(deltak(touch_idx>1)',p2,1); %later touches with corresponding tuning curve
        %FR2_all=simspikes(deltak(touch_idx>1)',p_all,1); %later touches with mixed tuning curve
        
        
        %ratioFR(av_unit,:)=[mean(FR(touch_idx==1)'./lambda1) 1/nanmean(FR(touch_idx==1)'./FR1_all(:)')];
        % response to first touch is X% higher than predicted
        %ratioFRb(av_unit,:)=nanmean(FR(touch_idx==1)'./FR1_all(:)');
        
        %% test if FR from first touch is higher than prediction
        [~,pvalFR(f)] = ttest(FR(touch_idx==1)',FR1_all(:)','Tail','right');
        %% test if FR is higher for the first touch
        [~,pval_atte(f)]=ttest2(FR(touch_idx==1),FR(touch_idx>1));
        
        FR_correct_incorrect_1(f,:)=[mean(FR_correct(touch_idx_correct==1)) mean(FR_incorrect(touch_idx_incorrect==1))];
        FR_correct_incorrect_later(f,:)=[mean(FR_correct(touch_idx_correct>1)) mean(FR_incorrect(touch_idx_incorrect>1))];
        
        % all trials        
        FR_all.correct_1=[FR_all.correct_1;FR_correct(touch_idx_correct==1)];
        FR_all.incorrect_1=[FR_all.incorrect_1;FR_incorrect(touch_idx_incorrect==1)];
        FR_all.correct_later=[FR_all.correct_later;FR_correct(touch_idx_correct>1)];
        FR_all.incorrect_later=[FR_all.incorrect_later;FR_incorrect(touch_idx_incorrect>1)];
        
        change_explained_by_adap(f)=min([(median(FR(touch_idx==1)-FR1_all(:)))/(median(FR(touch_idx==1))-median(FR(touch_idx>1))) 1e6]);
        
        adapt_mouse{this_mouse}= [adapt_mouse{this_mouse};change_explained_by_adap(f)];
        
        FR_touch(f,:)=[mean(FR(touch_idx==1)) mean(FR(touch_idx==2)) mean(FR(touch_idx==3)) mean(FR(touch_idx>3))]/mean(FR(touch_idx==1));
        FR_touch_correct(f,:)=[mean(FR_correct(touch_idx_correct==1)) mean(FR_correct(touch_idx_correct==2)) mean(FR_correct(touch_idx_correct==3)) mean(FR_correct(touch_idx_correct>3))]/mean(FR_correct(touch_idx_correct==1));
        FR_touch_incorrect(f,:)=[mean(FR_incorrect(touch_idx_incorrect==1)) mean(FR_incorrect(touch_idx_incorrect==2)) mean(FR_incorrect(touch_idx_incorrect==3)) mean(FR_incorrect(touch_idx_incorrect>3))]/mean(FR_incorrect(touch_idx_incorrect==1));
        
        FR_mouse{this_mouse}= [FR_mouse{this_mouse};FR_touch(f,:)];
        
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
errorbar(0.8,mean(FR_all.correct_1),std(FR_all.correct_1) ,'-.m')
errorbar(1.2,mean(FR_all.correct_later),std(FR_all.correct_later),'-.b')

% miss
errorbar(1.8,mean(FR_all.incorrect_1),std(FR_all.incorrect_1),'-.m')
errorbar(2.2,mean(FR_all.incorrect_later), std(FR_all.incorrect_later),'-.b')

ylabel('Firing rate [spikes/touch]')
box off
xlim([0.6 2.4])
xticks([1 2])
xticklabels({'Hit','Miss'})
plot([0.8 1.8],FR_correct_incorrect_1,'.r')
plot([1.2 2.2],FR_correct_incorrect_later,'.k')


%% ANOVA Fig 2C (FRs)
%% unbalanced
correct=[ones(size(FR_all.correct_1));zeros(size(FR_all.incorrect_1));ones(size(FR_all.correct_later));zeros(size(FR_all.incorrect_later))];
first_later=[ones(size(FR_all.correct_1));ones(size(FR_all.incorrect_1));ones(size(FR_all.correct_later))*2;ones(size(FR_all.incorrect_later))*2];
y=[FR_all.correct_1;FR_all.incorrect_1;FR_all.correct_later;FR_all.incorrect_later];


p = anovan(y,{correct,first_later},'display','off');

disp('-----------------------------------')
disp('Unbalanced anova First vs later responses (FR)')
disp(['p-value Effect hit vs miss = ' num2str(p(1))])
disp(['p-value Effect first vs later = ' num2str(p(2))])

%% balanced
nsample=size(FR_all.incorrect_1,1);
balanced_anova(FR_all.correct_1,FR_all.incorrect_1,FR_all.correct_later,FR_all.incorrect_later,nsample)


%% fraction between Responses in miss trials and hit trials.
disp(' ')
disp('-----------------------------------')
miss_hit_ratio=mean([FR_all.incorrect_1;FR_all.incorrect_later])/mean([FR_all.correct_1;FR_all.correct_later]);
disp(['Ration responses miss/hit = ' num2str(miss_hit_ratio)])

disp('-----------------------------------')

%% Change due to adaptation
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
