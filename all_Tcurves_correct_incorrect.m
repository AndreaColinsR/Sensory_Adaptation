function all_Tcurves_correct_incorrect
k3d=2;%1= k3d, 0=k_hor, 2=newk
counter=0;
av_unit=0;
%fig2=figure;
folder=[{'Glu32_19092017H'};{'Glu32_21092017H'};{'Glu43_22122017H'};{'Glu35_10112017H'};{'Glu35_13112017H_1'};{'Glu35_13112017H_2'}];
%folder2=[{'Glu32_190917a_20170919_'},{'Glu32_210917a_20170921_'},{'Glu43_221217a_20171222_'},{'Glu35_101117a_20171110_'},{'Glu35_131117a_20171113_'},{'Glu35_131117a_20171113_'}];
[~,txt]=xlsread('Adaptation_units_list.xlsx',1);
%colours=[{'b'};{'g'};{'r'};{'c'};{'m'};{'k'}];
deltak_correct_all_trials=[];
deltak_incorrect_all_trials=[];
deltak_correct_later=[];
deltak_incorrect_later=[];
for f=1:size(folder,1)
    
    cd(folder{f})
    info=xlsread('ledtrials.xlsx');
    go_trials=find(info(:,7)==2);
    idx_correct=info(go_trials,6)==1;
    idx_incorrect=info(go_trials,6)==2;
    %cd('./modelling/')
    
    if k3d==1
        %k3D
        %         load('deltak_c2.mat')
        %         k_c2=touches_matrix/0.047;
        %         load('deltak_c1.mat')
        %         k_c1=touches_matrix/0.047;
        [~,~,~,matrix_norm]=normalised_by_trial(15,0);
        k_c1=matrix_norm(:,:,1);
        k_c2=matrix_norm(:,:,2);
    elseif k3d==0
        %k_h
        load('kinematicsgoc2.mat')
        k_c2=azimuth;
        load('kinematicsgoc1.mat')
        k_c1=azimuth;
        
        %         for i=1:size(azimuth,1)
        %             k_c1(i,:)=[0 abs(diff(k_c1(i,:)))];
        %            k_c2(i,:)=[0 abs(diff(k_c2(i,:)))];
        %         end
    else
        %newk
        %          load('surprise.mat')
        %          k_c2=[zeros(size(surprise,1),500) surprise];
        %          k_c1=[zeros(size(surprise,1),500) surprise];
        
        load('newk1.mat')
        k_c1=newk1*1.7/0.047;
        load('newk2.mat')
        k_c2=newk2/0.047;
        
        %           k_c1=distance_matrix(:,:,1)*0.047;
        %           k_c2=distance_matrix(:,:,2)*0.047;
        
        %         load('kinematicsgoc2.mat')
        %         k_c2=k_ver/0.047;
        %         load('kinematicsgoc1.mat')
        %         k_c1=k_ver*1.7/0.047;
        
        %              load('roll_c1.mat')
        %              k_c1=touches_matrix;
        %               load('roll_c2.mat')
        %             k_c2=touches_matrix;
        
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
    
    for unit=1:size(ff,1)
        onthelist(txt,[folder{f} ff(unit).name])
        
        if onthelist(txt,[folder{f} ff(unit).name])
            load(ff(unit).name)
            total_psth=psth;
            total_psth_correct=psth(idx_correct,:);
            total_psth_incorrect=psth(idx_incorrect,:);
            delete=size(total_psth,2)-size(k_c1,2);
            k_c1=[zeros(size(k_c1,1),delete),k_c1];
            k_c2=[zeros(size(k_c1,1),delete),k_c2];
            clear deltak touch_idx FR
            
            [x1,p1,x2,p2,x1_norm,p1_norm,x2_norm,p2_norm,x_all,p_all,deltak,touch_idx,FR,~,FR_prev,whisker]=Tcurve_touches_dk(touches_matrix,touches_whisker,total_psth,k_c1,k_c2,0,4);
            [x1_correct,p1_correct,x2_correct,p2_correct,~,~,~,~,~,~,deltak_correct,touch_idx_correct,FR_correct,~,~,~]=Tcurve_touches_dk(touches_matrix_correct,touches_whisker_correct,total_psth_correct,k_c1_correct,k_c2_correct,0,4);
            [x1_incorrect,p1_incorrect,x2_incorrect,p2_incorrect,~,~,~,~,~,~,deltak_incorrect,touch_idx_incorrect,FR_incorrect,~,~,~]=Tcurve_touches_dk(touches_matrix_incorrect,touches_whisker_incorrect,total_psth_incorrect,k_c1_incorrect,k_c2_incorrect,0,4);
            
            [lambda1,~]=Tcurve_touches_dk_cross_val(touches_matrix,touches_whisker,total_psth,k_c1,k_c2,0,4);
            if ~isnan(p1)
                av_unit=av_unit+1;
                
                spikes_incorrect_1(av_unit)=sum(FR_incorrect(touch_idx_incorrect==1));
                coeff1(av_unit,:)=p1;
                coeff2(av_unit,:)=p2;
                %figure(fig2)
                f1 = polyval(p1,x1);
                f2 = polyval(p2,x1);
                f1_correct = polyval(p1_correct,x1_correct);
                f2_correct = polyval(p2_correct,x1_correct);
                f1_incorrect = polyval(p1_incorrect,x1_incorrect);
                f2_incorrect = polyval(p2_incorrect,x1_incorrect);
                subplot(3,3,1)
                f_pred=sqrt(f1*FR_prev(1));
                %plot([-0.01; x1],[FR_prev(1); f1],'-+g')
                hold on
                %plot([-0.01; x1],[FR_prev(2); f2],'-ok')
                plot(x1_correct,f1_correct,'-+b')
                plot(x1_incorrect,f1_incorrect,'-+r')
                f_all = polyval(p_all,x1);
                %pause
                % plot(x1,f_all,'-oc')
                %R=coeff_det(f1/p1(2),f2/p2(2));
                pos_slope(av_unit)=p_all(1);
                
                subplot(3,3,2)
                %plot(p2_correct(1),p1_correct(1),'ob')
                hold on
                %plot(p2_incorrect(1),p1_incorrect(1),'or')
                %plot([p2_correct(1) p2_incorrect(1)],[p1_correct(1)
                %p1_incorrect(1)],'k')
                %plot([1 2],[p1_correct(1) p1_incorrect(1)],'r')
                %plot([1 2],[p2_correct(1) p2_incorrect(1)],'k')
                %plot(p2(1),p_all(1),'or')
                
                
                
                subplot(3,3,3)
                %plot(p2_correct(2),p1_correct(2),'ob')
                hold on
                %plot(p2_incorrect(2),p1_incorrect(2),'or')
                %plot([p2_correct(2) p2_incorrect(2)],[p1_correct(2) p1_incorrect(2)],'k')
                %plot([1 2],[p2_correct(2), p2_incorrect(2)],'k')
                %plot([1 2],[p1_correct(2), p1_incorrect(2)],'r')
                subplot(3,3,4)
                %FR1_1=simspikes(deltak(touch_idx==1)',p1,1); %first touch with correspondent tuning curve
                FR1_all=simspikes(deltak(touch_idx==1)',p2,1); %first touch with mixed tunning cuve
                FR2=simspikes(deltak(touch_idx>1)',p2,1); %later touches with corresponding tuning curve
                FR2_all=simspikes(deltak(touch_idx>1)',p_all,1); %later touches with mixed tuning curve
                
                hold on
                xk(av_unit,:)=[mean(FR(touch_idx==1)),mean(FR1_all(:))];
                xr(av_unit,:)=[mean(FR(touch_idx==1)),mean(lambda1(:))];
                plot(mean(FR(touch_idx==1)),mean(lambda1(:)),'or')
                plot(mean(FR(touch_idx==1)),mean(FR1_all(:)),'ok')
                %plot(mean(abs(FR(touch_idx==1)-lambda1')),mean(abs(FR(touch_idx==1)-FR1_all)),'ok')
                
                subplot(3,3,5)
                
                %errorbar([1 2 3 4],[mean(FR1_1(:)) mean(FR1_all(:)) mean(FR2_all(:)) mean(FR2(:))],[std(FR1_1(:)) std(FR1_all(:)) std(FR2_all(:)) std(FR2(:))],'.-')
                plot([1 2],[mean(FR(touch_idx==1)'./lambda1) 1/nanmean(FR(touch_idx==1)'./FR1_all(:)')],'bo-')
                ratioFR(av_unit,:)=[mean(FR(touch_idx==1)'./lambda1) 1/nanmean(FR(touch_idx==1)'./FR1_all(:)')];
                % response to first touch is X% higher than predicted
                ratioFRb(av_unit,:)=nanmean(FR(touch_idx==1)'./FR1_all(:)');
                hold on
                
                %% test if FR from first touch is higher than prediction
                [hFR(av_unit),pvalFR(av_unit)] = ttest(FR(touch_idx==1)',FR1_all(:)','Tail','right');
                
                
                %subplot(3,3,7)
                %plot([1 2],[1 mean(FR2(:))/mean(FR2_all(:))],'bo-')
                
                error1(av_unit)=mean(abs(FR(touch_idx==1)-lambda1(:)));
                error2(av_unit)=mean(abs(FR(touch_idx==1)-FR1_all(:)));
                %bar(av_unit,error1,'r')
                %hold on
                %bar(av_unit,-error2,'k')
                ratioFR2(av_unit,:)=[1 mean(FR2(:))/mean(FR2_all(:))];
                
                %
                subplot(3,3,8)
                %plot([ 1 2],[FR_prev(1) mean(FR(touch_idx==1))],'r-o')
                %plot([ 1 2],[FR_prev(2) mean(FR(touch_idx>1))],'k-o')
                hold on
                FR_correct_incorrect_1(av_unit,:)=[mean(FR_correct(touch_idx_correct==1)) mean(FR_incorrect(touch_idx_incorrect==1))];
                FR_correct_incorrect_later(av_unit,:)=[mean(FR_correct(touch_idx_correct>1)) mean(FR_incorrect(touch_idx_incorrect>1))];
                
                plot([1 2],FR_correct_incorrect_1(av_unit,:),'-r')
                plot([1 2],FR_correct_incorrect_later(av_unit,:),'-k')
               title('FR')
                %deltak_correct,touch_idx_correct,FR_correct
                
                %deltak_touch(f,:)=[mean(deltak(touch_idx==1)) mean(delta2) mean(delta3) mean(delta4)];
                %plot([1 2 3 4],deltak_touch(f,:))
                
                
                subplot(3,3,9)
                %hold on
                %errorbar([1 2],[mean(FR((touch_idx>1) & (whisker==1))),mean(FR((touch_idx>1) & (whisker==2)))],[std(FR((touch_idx>1) & (whisker==1))),std(FR((touch_idx>1) & (whisker==2)))])
                hw(av_unit)=ttest2(FR((whisker==1)),FR((whisker==2)));
                
                %         plot(FR_prev(1),p1(2),['or'])
                %         hold on
                %         plot(FR_prev(2),p2(2),['dk'])
                %         xlabel('preactivity')
                %         ylabel('Intercept')
                
                %          1/nanmean(FR(touch_idx==1)'./FR1_all(:)')
                %         nanmean(FR(touch_idx==1)-FR1_all(:))
                %mean(FR(touch_idx==2))
                %(mean(FR(touch_idx==1))-mean(FR(touch_idx==2)))
                %change_explained_by_adap(av_unit)=(1-ratioFR(av_unit,2))/(1-(mean(FR(touch_idx>1))/mean(FR(touch_idx==1))))
                [h_atte(av_unit),pval_atte(av_unit)]=ttest2(FR(touch_idx==1),FR(touch_idx>1));
                
                
                if h_atte(av_unit) %(mean(FR(touch_idx==1))-mean(FR(touch_idx>1)))>0
                    change_explained_by_adap(av_unit)=min([(median(FR(touch_idx==1)-FR1_all(:)))/(mean(FR(touch_idx==1))-mean(FR(touch_idx>1))) 1e6]);
                    
                    %change_explained_by_adap(av_unit)=(1-ratioFR(av_unit,2))/(1-(mean(FR(touch_idx>3))/mean(FR(touch_idx==1))));
                    
                else
                    change_explained_by_adap(av_unit)=NaN;
                end
                FR_touch(av_unit,:)=[mean(FR(touch_idx==1)) mean(FR(touch_idx==2)) mean(FR(touch_idx==3)) mean(FR(touch_idx>3))]/mean(FR(touch_idx==1));
                FR_touch_correct(av_unit,:)=[mean(FR_correct(touch_idx_correct==1)) mean(FR_correct(touch_idx_correct==2)) mean(FR_correct(touch_idx_correct==3)) mean(FR_correct(touch_idx_correct>3))]/mean(FR_correct(touch_idx_correct==1));
                FR_touch_incorrect(av_unit,:)=[mean(FR_incorrect(touch_idx_incorrect==1)) mean(FR_incorrect(touch_idx_incorrect==2)) mean(FR_incorrect(touch_idx_incorrect==3)) mean(FR_incorrect(touch_idx_incorrect>3))]/mean(FR_incorrect(touch_idx_incorrect==1));
                
                %                 plot([1 2 3 4],FR_touch(av_unit,:),'k')
                %                 hold on
                %                 plot([1 2 3 4],FR_touch_correct(av_unit,:),'b')
                %                 plot([1 2 3 4],FR_touch_incorrect(av_unit,:),'r')
                
                
                
                if p2(1)<p1(1)
                    counter=counter+1;
                end
                
            end
            
            clear FR FR_prev FR1_all FR2_all FR2
        end
    end
    
    subplot(3,3,6)
    deltak_touch(f,:)=[mean(deltak(touch_idx==1)) mean(deltak(touch_idx==2)) mean(deltak(touch_idx==3)) mean(deltak(touch_idx>3))];
    deltak_touch_correct(f,:)=[mean(deltak_correct(touch_idx_correct==1)) mean(deltak_correct(touch_idx_correct==2)) mean(deltak_correct(touch_idx_correct==3)) mean(deltak_correct(touch_idx_correct>3))];
    deltak_touch_incorrect(f,:)=[mean(deltak_incorrect(touch_idx_incorrect==1)) mean(deltak_incorrect(touch_idx_incorrect==2)) mean(deltak_incorrect(touch_idx_incorrect==3)) mean(deltak_incorrect(touch_idx_incorrect>3))];
    
    hold on
    
    subplot(3,3,7)
    deltak_correct_incorrect_1(f,:)=[mean(deltak_correct(touch_idx_correct==1)) mean(deltak_incorrect(touch_idx_incorrect==1))];
    deltak_correct_incorrect_later(f,:)=[mean(deltak_correct(touch_idx_correct>1)) mean(deltak_incorrect(touch_idx_incorrect>1))];
    
    deltak_correct_all_trials=[deltak_correct_all_trials deltak_correct(touch_idx_correct==1)];
    deltak_incorrect_all_trials=[deltak_incorrect_all_trials deltak_incorrect(touch_idx_incorrect==1)];
    
    deltak_correct_later=[deltak_correct_later deltak_correct(touch_idx_correct>1)];
    deltak_incorrect_later=[deltak_incorrect_later deltak_incorrect(touch_idx_incorrect>1)];
    
    
    plot([1 2],deltak_correct_incorrect_1,'-r')
    hold on
    plot([1 2],deltak_correct_incorrect_later,'-k')
    title('delta k')
    
    clear deltak touch_idx
    cd ..

end

subplot(3,3,1)
box off
xlabel('\Delta \kappa')
ylabel('Normalised FR')
title('Tuning curves')


subplot(3,3,2)
hold on
% plot([-25 600],[-25 600],'-b','Linewidth',2)
% plot([0 0],[-25 600],'-k','Linewidth',1)
% plot([-25 600],[0 0],'-k','Linewidth',1)
% box off
% xlabel('Later touches slope')
% ylabel('First touches slope')
% ylim([-25 600])
% xlim([-25 600])

% plot(mean(coeff2(:,1)),mean(coeff1(:,1)),'.k')
% errorbar(mean(coeff2(:,1)),mean(coeff1(:,1)),std(coeff2(:,1)),'horizontal','.k')
% errorbar(mean(coeff2(:,1)),mean(coeff1(:,1)),std(coeff1(:,1)),'vertical','.k')
subplot(3,3,3)
box off
% xlabel('Later touches intercept')
% ylabel('First touches intercept')
% ylim([0 6])
% xlim([0 6])
% plot([0 6],[0 6],'-b','Linewidth',2)
% plot(mean(coeff2(:,2)),mean(coeff1(:,2)),'.k')
% errorbar(mean(coeff2(:,2)),mean(coeff1(:,2)),std(coeff2(:,2)),'horizontal','.k')
% errorbar(mean(coeff2(:,2)),mean(coeff1(:,2)),std(coeff1(:,2)),'vertical','.k')
%Does the intercept significantly decrease across the population
[hintercept,pvalintercept] = ttest(coeff1(:,2)',coeff2(:,2)','Tail','right');
% how much decrease
intercept_decrease=1-mean(coeff2(:,2)./coeff1(:,2));
%Does the slope significantly decrease across the population
[hslope,pvalslope] = ttest(coeff1(:,1),coeff2(:,1),'Tail','right');
slope_decrease=coeff2(:,1)./coeff1(:,1);

subplot(3,3,4)
box off
xlabel('Mean FR')
ylabel('Predicted mean FR')
plot([0 8],[0 8],'b')
xlim([0 8])
ylim([0 8])
%legend('First touch TC','Later touches TC')
title('Prediction of first touch responses')
plot(mean(xk(:,1)),mean(xk(:,2)),'.k')
plot(mean(xr(:,1)),mean(xr(:,2)),'.r')
errorbar(mean(xk(:,1)),mean(xk(:,2)),std(xk(:,1)),'horizontal','.k')
errorbar(mean(xk(:,1)),mean(xk(:,2)),std(xk(:,2)),'vertical','.k')
errorbar(mean(xr(:,1)),mean(xr(:,2)),std(xr(:,1)),'horizontal','.r')
errorbar(mean(xr(:,1)),mean(xr(:,2)),std(xr(:,2)),'vertical','.r')


subplot(3,3,5)
box off
xlim([0.8 2.2])
ylim([0 2])
plot([1 2],[1 1],'k')
xticks([1 2])
ylabel('Ratio FR tuning curves')
xticklabels({'later touch ratio','later touches ratio'})
errorbar([1 2], nanmean(ratioFR,1),nanstd(ratioFR,[],1),'.-k','LineWidth',2)


subplot(3,3,6)
xlabel('Touch number')
xticks([1 2 3 4])
ylabel('K')
xlim([0.5 4.5])

errorbar(1:4,mean(deltak_touch,1),std(deltak_touch,[],1),'k')
errorbar(1:4,mean(deltak_touch_correct,1),std(deltak_touch_correct,[],1),'b')
errorbar(1:4,mean(deltak_touch_incorrect,1),std(deltak_touch_incorrect,[],1),'r')

[hd1,pd1] = ttest2(deltak_touch(:,1),deltak_touch(:,2));
[hd2,pd2] = ttest2(deltak_touch(:,2),deltak_touch(:,3));
[hd3,pd3] = ttest2(deltak_touch(:,3),deltak_touch(:,4));

subplot(3,3,8)
box off
xlim([0.8 2.2])
errorbar([1 2],mean(FR_correct_incorrect_1),std(FR_correct_incorrect_1),'-m')
errorbar([1 2],mean(FR_correct_incorrect_later),std(FR_correct_incorrect_later),'-b')
%%% correct vs incorrect
[h_FR_1,pFR_1] = ttest2(FR_correct_incorrect_1(:,1),FR_correct_incorrect_1(:,2));
[h_FR_later,pFR_later] = ttest2(FR_correct_incorrect_later(:,1),FR_correct_incorrect_later(:,2));

%%% 1st vs later
[h_FR_correct,pFR_correct] = ttest2(FR_correct_incorrect_1(:,1),FR_correct_incorrect_later(:,1));
[h_FR_incorrect,pFR_incorrect] = ttest2(FR_correct_incorrect_1(:,2),FR_correct_incorrect_later(:,2));


subplot(3,3,7)
box off
xlim([0.8 2.2])
errorbar([1 2],mean(deltak_correct_incorrect_1),std(deltak_correct_incorrect_1),'-m')
errorbar([1 2],mean(deltak_correct_incorrect_later),std(deltak_correct_incorrect_later),'-b')

%%% correct vs incorrect
[h_dk_1,pdk_1] = ttest2(deltak_correct_incorrect_1(:,1),deltak_correct_incorrect_1(:,2));
[h_dk_later,pdk_later] = ttest2(deltak_correct_incorrect_later(:,1),deltak_correct_incorrect_later(:,2));

[h_dk_1,pdk_1] = ttest2(deltak_correct_all_trials,deltak_incorrect_all_trials)
[h_dk_later,pdk_later] = ttest2(deltak_correct_later,deltak_incorrect_later)
size(deltak_incorrect_all_trials)

%%% 1st vs later
[h_dk_correct,pdk_correct] = ttest2(deltak_correct_incorrect_1(:,1),deltak_correct_incorrect_later(:,1));
[h_dk_incorrect,pdk_incorrect] = ttest2(deltak_correct_incorrect_1(:,2),deltak_correct_incorrect_later(:,2));

[h_dk_correct,pdk_correct] = ttest2(deltak_correct_all_trials,deltak_correct_later)
[h_dk_incorrect,pdk_incorrect] = ttest2(deltak_incorrect_all_trials,deltak_incorrect_later)

pause

subplot(3,3,9)
hold on
errorbar(1:4,mean(FR_touch,1),std(FR_touch,[],1),'k')
errorbar(1:4,mean(FR_touch_correct,1),std(FR_touch_correct,[],1),'b')
errorbar(1:4,mean(FR_touch_incorrect,1),std(FR_touch_incorrect,[],1),'r')


[h1,p1] = ttest(FR_touch(:,2),1,'Tail','left');
[h2,p2] = ttest2(FR_touch(:,2),FR_touch(:,3),'Tail','right');
[h3,p3] = ttest2(FR_touch(:,3),FR_touch(:,4),'Tail','right');

% how many units show adaptation?
lower=sum(hFR);

% change_explained_by_adap;
% nanmean(change_explained_by_adap)
% nanstd(change_explained_by_adap)


%subplot(3,3,9)
% FO=fit(X',Y','poly1');
% plot([0 max(X)], FO.p1*[0 max(X)]+FO.p2,'k')
% response to first touch is X% higher than predicted
% mean(ratioFRb(ratioFRb<10))
% mean(pval_atte(h_atte>0))
clear deltak touch_idx
 nanmedian(spikes_incorrect_1)
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