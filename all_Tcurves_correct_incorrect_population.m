function [FR_touch,FR_touch_correct,FR_touch_incorrect]=all_Tcurves_correct_incorrect_population
k3d=2;%1= k3d, 0=k_hor, 2=newk
counter=0;
av_unit=0;
folder=[{'Glu32_19092017H'};{'Glu32_21092017H'};{'Glu43_22122017H'};{'Glu35_10112017H'};{'Glu35_13112017H_1'};{'Glu35_13112017H_1'}];
[~,txt]=xlsread('Adaptation_units_list.xlsx',1);
percentile_1=0.95;
percentile_2=0.95;

FR_correct_1_all=[];
FR_incorrect_1_all=[];
FR_correct_later_all=[];
FR_incorrect_later_all=[];

for f=1:size(folder,1)
    
    cd(folder{f})
    info=xlsread('ledtrials.xlsx');
    go_trials=find(info(:,7)==2);
    idx_correct=info(go_trials,6)==1;
    idx_incorrect=info(go_trials,6)==2;
    
    if k3d==1
        [~,~,~,matrix_norm]=normalised_by_trial(15,0);
        k_c1=matrix_norm(:,:,1);
        k_c2=matrix_norm(:,:,2);
    elseif k3d==0
        %k_h
        load('kinematicsgoc2.mat')
        k_c2=azimuth;
        load('kinematicsgoc1.mat')
        k_c1=azimuth;
        
    else
        
        
        load('newk1.mat')
        k_c1=newk1*1.7/0.047;
        load('newk2.mat')
        k_c2=newk2/0.047;
        
        
    end
    
    
    load('all_touches.mat')
    load('touches_whisker.mat')
    ff=dir('*neural_data_S*');
    %% select correct and incorrect trials
    touches_matrix_correct=touches_matrix(idx_correct,:);
    touches_whisker_correct=touches_whisker(idx_correct,:,:);
    k_c1_correct=k_c1(idx_correct,:);
    k_c2_correct=k_c2(idx_correct,:);
    touches_matrix_incorrect=touches_matrix(idx_incorrect,:);
    touches_whisker_incorrect=touches_whisker(idx_incorrect,:,:);
    k_c1_incorrect=k_c1(idx_incorrect,:);
    k_c2_incorrect=k_c2(idx_incorrect,:);
    total_psth=zeros(size(k_c1));
    total_psth_correct=zeros(size(k_c1_correct));
    total_psth_incorrect=zeros(size(k_c1_incorrect));
    
    for unit=1:size(ff,1)
        
        if onthelist(txt,[folder{f} ff(unit).name])
            
            load(ff(unit).name)
            total_psth=psth+total_psth;
            total_psth_correct=psth(idx_correct,:)+total_psth_correct;
            total_psth_incorrect=psth(idx_incorrect,:)+total_psth_incorrect;
            delete=size(total_psth,2)-size(k_c1,2);
            k_c1=[zeros(size(k_c1,1),delete),k_c1];
            k_c2=[zeros(size(k_c1,1),delete),k_c2];
            clear deltak touch_idx FR
        end
    end
    
    [x1,p1,x2,p2,x1_norm,p1_norm,x2_norm,p2_norm,x_all,p_all,deltak,touch_idx,FR,~,FR_prev,whisker]=Tcurve_touches_dk(touches_matrix,touches_whisker,total_psth,k_c1,k_c2,0,4,0,percentile_1,percentile_2);
    [x1_correct,p1_correct,x2_correct,p2_correct,~,~,~,~,~,~,~,touch_idx_correct,FR_correct,~,~,~,spikes_correct]=Tcurve_touches_dk(touches_matrix_correct,touches_whisker_correct,total_psth_correct,k_c1_correct,k_c2_correct,0,4,0,percentile_1,percentile_2);
    [x1_incorrect,p1_incorrect,x2_incorrect,p2_incorrect,~,~,~,~,~,~,~,touch_idx_incorrect,FR_incorrect,~,~,~,spikes_incorrect]=Tcurve_touches_dk(touches_matrix_incorrect,touches_whisker_incorrect,total_psth_incorrect,k_c1_incorrect,k_c2_incorrect,0,4,0,percentile_1,percentile_2);
    
    if f==1
    raster_plot(spikes_correct,touch_idx_correct,[1 2 4 5])
    raster_plot(spikes_incorrect,touch_idx_incorrect,[7 8 10 11])
    end
    
    [lambda1,~]=Tcurve_touches_dk_cross_val(touches_matrix,touches_whisker,total_psth,k_c1,k_c2,0,4);
    
    if ~isnan(p1_incorrect)
        av_unit=av_unit+1;
        
        spikes_incorrect_1(av_unit)=sum(FR_correct(touch_idx_correct==1));
        coeff1(av_unit,:)=p1;
        coeff2(av_unit,:)=p2;
        
        
        f1_correct = polyval(p1_correct,x1_correct);
        f2_correct = polyval(p2_correct,x1_correct);
        f1_incorrect = polyval(p1_incorrect,x1_incorrect);
        f2_incorrect = polyval(p2_incorrect,x1_incorrect);
        pos_slope(av_unit)=p_all(1);
        
        
        slope_first(av_unit,:)=[p1_correct(1) p1_incorrect(1)];
        slope_later(av_unit,:)=[p2_correct(1) p2_incorrect(1)];
        intercept_first(av_unit,:)=[p1_correct(2), p1_incorrect(2)];
        intercept_later(av_unit,:)=[p2_correct(2), p2_incorrect(2)];
        
        
        FR1_all=simspikes(deltak(touch_idx==1)',p2,1); %first touch with mixed tunning cuve
        FR2=simspikes(deltak(touch_idx>1)',p2,1); %later touches with corresponding tuning curve
        FR2_all=simspikes(deltak(touch_idx>1)',p_all,1); %later touches with mixed tuning curve
        
        
        ratioFR(av_unit,:)=[mean(FR(touch_idx==1)'./lambda1) 1/nanmean(FR(touch_idx==1)'./FR1_all(:)')];
        % response to first touch is X% higher than predicted
        ratioFRb(av_unit,:)=nanmean(FR(touch_idx==1)'./FR1_all(:)');
        
        %% test if FR from first touch is higher than prediction
        [hFR(av_unit),pvalFR(av_unit)] = ttest(FR(touch_idx==1)',FR1_all(:)','Tail','right');
        
        
        error1(av_unit)=mean(abs(FR(touch_idx==1)-lambda1(:)));
        error2(av_unit)=mean(abs(FR(touch_idx==1)-FR1_all(:)));
        
        ratioFR2(av_unit,:)=[1 mean(FR2(:))/mean(FR2_all(:))];
        
        
        FR_correct_incorrect_1(av_unit,:)=[mean(FR_correct(touch_idx_correct==1)) mean(FR_incorrect(touch_idx_incorrect==1))];
        FR_correct_incorrect_later(av_unit,:)=[mean(FR_correct(touch_idx_correct>1)) mean(FR_incorrect(touch_idx_incorrect>1))];
        
        %all trials
        FR_correct_1_all=[FR_correct_1_all;FR_correct(touch_idx_correct==1)];
        FR_incorrect_1_all=[FR_incorrect_1_all;FR_incorrect(touch_idx_incorrect==1)];
        FR_correct_later_all=[FR_correct_later_all;FR_correct(touch_idx_correct>1)];
        FR_incorrect_later_all=[FR_incorrect_later_all;FR_incorrect(touch_idx_incorrect>1)];
        
        
        hw(av_unit)=ttest2(FR((whisker==1)),FR((whisker==2)));
        
        [h_atte(av_unit),pval_atte(av_unit)]=ttest2(FR(touch_idx==1),FR(touch_idx>1));
        
        
        if h_atte(av_unit) %(mean(FR(touch_idx==1))-mean(FR(touch_idx>1)))>0
            change_explained_by_adap(av_unit)=min([(median(FR(touch_idx==1)-FR1_all(:)))/(median(FR(touch_idx==1))-median(FR(touch_idx>1))) 1e6]);
            
            %change_explained_by_adap(av_unit)=(1-ratioFR(av_unit,2))/(1-(mean(FR(touch_idx>3))/mean(FR(touch_idx==1))));
            
        else
            change_explained_by_adap(av_unit)=NaN;
        end
        
        FR_touch(av_unit,:)=[mean(FR(touch_idx==1)) mean(FR(touch_idx==2)) mean(FR(touch_idx==3)) mean(FR(touch_idx>3))]/mean(FR(touch_idx==1));
        FR_touch_correct(av_unit,:)=[mean(FR_correct(touch_idx_correct==1)) mean(FR_correct(touch_idx_correct==2)) mean(FR_correct(touch_idx_correct==3)) mean(FR_correct(touch_idx_correct>3))]/mean(FR_correct(touch_idx_correct==1));
        FR_touch_incorrect(av_unit,:)=[mean(FR_incorrect(touch_idx_incorrect==1)) mean(FR_incorrect(touch_idx_incorrect==2)) mean(FR_incorrect(touch_idx_incorrect==3)) mean(FR_incorrect(touch_idx_incorrect>3))]/mean(FR_incorrect(touch_idx_incorrect==1));
        
        
        if p2(1)<p1(1)
            counter=counter+1;
        end
        
        clear FR FR_prev FR1_all FR2_all FR2
    end
    
    clear touch_idx
    cd ..
    
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


[h_FR_1,pFR_1] = ttest2(FR_correct_1_all,(FR_incorrect_1_all));
[h_FR_later,pFR_later] = ttest2(FR_correct_later_all,FR_incorrect_later_all);

%
[h_FR_correct,pFR_correct] = ttest2(FR_correct_1_all,FR_correct_later_all);
[h_FR_incorrect,pFR_incorrect] = ttest2(FR_incorrect_1_all,FR_incorrect_later_all);
%pause

%%ANOVA Fig 2C (FRs)
%% unbalanced
correct=[ones(size(FR_correct_1_all));zeros(size(FR_incorrect_1_all));ones(size(FR_correct_later_all));zeros(size(FR_incorrect_later_all))];
first_later=[ones(size(FR_correct_1_all));ones(size(FR_incorrect_1_all));ones(size(FR_correct_later_all))*2;ones(size(FR_incorrect_later_all))*2];
y=[FR_correct_1_all;FR_incorrect_1_all;FR_correct_later_all;FR_incorrect_later_all];
p = anovan(y,{correct,first_later});
p(1)
p(2)

mean([FR_incorrect_1_all;FR_incorrect_later_all])/mean([FR_correct_1_all;FR_correct_later_all])
%pause

%%
%balanced
% nsample=size(FR_incorrect_1_all,1);
% correct=[ones(nsample,1);zeros(nsample,1);ones(nsample,1);zeros(nsample,1)];
% first_later=[ones(nsample,1);ones(nsample,1);ones(nsample,1)*2;ones(nsample,1)*2];
% for i=1:100
% y=[FR_correct_1_all(randperm(numel(FR_correct_1_all),nsample));...
%     FR_incorrect_1_all(randperm(numel(FR_incorrect_1_all),nsample));...
%     FR_correct_later_all((randperm(numel(FR_correct_later_all),nsample)));...
%     FR_incorrect_later_all((randperm(numel(FR_incorrect_later_all),nsample)))];
% p(:,i) = anovan(y,{correct,first_later},'display','off');
% clear y
% end
%
%  [mean(p(1,:)) std(p(1,:))]
%  [mean(p(2,:)) std(p(2,:))]

%%%%% ttest FR
%errorbar([1 2],mean(FR_correct_incorrect_1),std(FR_correct_incorrect_1),'-m')
%errorbar([1 2],mean(FR_correct_incorrect_later),std(FR_correct_incorrect_later),'-b')
%[h_FR_1,pFR_1] = ttest2(FR_correct_incorrect_1(:,1),FR_correct_incorrect_1(:,2));
%[h_FR_later,pFR_later] = ttest2(FR_correct_incorrect_later(:,1),FR_correct_incorrect_later(:,2));



[h1,p1] = ttest(FR_touch(:,2),1,'Tail','left');
[h2,p2] = ttest2(FR_touch(:,2),FR_touch(:,3),'Tail','right');
[h3,p3] = ttest2(FR_touch(:,3),FR_touch(:,4),'Tail','right');

% how many units show adaptation?
lower=sum(hFR);
%mean(pvalFR(hFR>0))

%subplot(3,3,9)
% FO=fit(X',Y','poly1');
% plot([0 max(X)], FO.p1*[0 max(X)]+FO.p2,'k')
% response to first touch is X% higher than predicted
% mean(ratioFRb(ratioFRb<10))
% mean(pval_atte(h_atte>0))
% clear deltak touch_idx
% min(spikes_incorrect_1)

% figure
% histogram(change_explained_by_adap,0:0.1:1.5)
change_explained_by_adap
mean(change_explained_by_adap)
std(change_explained_by_adap)
end

function FR=simspikes(x,p,Nsim)
lambda=p(1)*x+p(2);
lambda=repmat(lambda,1,Nsim);
FR = lambda;
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