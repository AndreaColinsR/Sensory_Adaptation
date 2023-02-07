function [change_explained_by_adap,FR_touch]=all_Tcurves(fig4)
k3d=2;%1= k3d, 0=k_hor, 2=newk
do_extra_plot=0;
counter=0;
counter2=1;
av_unit=0;
percentile_1=0.95;
percentile_2=0.95;
NpointsTcurve=4;

folder=[{'Glu32_19092017H'};{'Glu32_21092017H'};{'Glu43_22122017H'};{'Glu35_10112017H'};{'Glu35_13112017H_1'};{'Glu35_13112017H_2'}];
[~,txt]=xlsread('Adaptation_units_list.xlsx',1);

nunits=size(txt,1);
FR_touch=nan(nunits,4);
for f=1:size(folder,1)
    
    cd(folder{f})
    
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
    clear total_psth
    total_psth=zeros(size(touches_matrix));
    for unit=1:size(ff,1)
        
        if onthelist(txt,[folder{f} ff(unit).name])
            load(ff(unit).name)
            total_psth=psth;
            
            
            delete=size(total_psth,2)-size(k_c1,2);
            k_c1=[zeros(size(k_c1,1),delete),k_c1];
            k_c2=[zeros(size(k_c1,1),delete),k_c2];
            %     total_psth(:,1:delete)=[];
            %     touches_matrix(:,1:delete)=[];
            %     touches_whisker(:,1:delete)=[];
            %     size(total_psth)
            %     size(k_c1)
            clear deltak touch_idx FR
            
            [x1,p1,x2,p2,x1_norm,p1_norm,x2_norm,p2_norm,x_all,p_all,deltak,touch_idx,FR,~,FR_prev,whisker]=Tcurve_touches_dk(touches_matrix,touches_whisker,total_psth,k_c1,k_c2,do_extra_plot,NpointsTcurve,0,percentile_1,percentile_2);
            deltak=abs(deltak);
            
            % [x1,p1,x2,p2,x1_norm,p1_norm,x2_norm,p2_norm,x_all,p_all,deltak,touch_idx,FR,deltak_norm,y2,y2_norm]=Tcurve_touches_dk_norm_bywhisker(touches_matrix,touches_whisker,total_psth,k_c1,k_c2,1,6);
            %[~,~,~,~,~,~,~,~,x_all,p_all,~,~,~,~,~,~]=Tcurve_slips_dk(slips_matrix,slips_whisker,total_psth,k_c1,k_c2,0,4);
            [lambda1,~]=Tcurve_touches_dk_cross_val(touches_matrix,touches_whisker,total_psth,k_c1,k_c2,0,4);
            
            if ~isnan(x1)
                av_unit=av_unit+1;
                
                coeff1(av_unit,:)=p1;
                coeff2(av_unit,:)=p2;
                
                f1 = polyval(p1,x1);
                f2 = polyval(p2,x1);
                figure(fig4)
                if f==4 && counter2==1
                    
                    subplot(2,4,2)
                    nsamples2=plotTCurve(deltak(touch_idx>1),FR(touch_idx>1),'k',NpointsTcurve,percentile_2);
                    nsamples1=plotTCurve(deltak(touch_idx==1),FR(touch_idx==1),'r',NpointsTcurve,percentile_1);
                    xlabel('|\Delta \kappa| [1/mm]')
                    ylabel('FR [spikes/touch]')
                    box off
                    counter2=2;
                end
                
                subplot(2,4,3)
                plot([-0.01; x1],[FR_prev(1); f1],'-or')
                hold on
                plot([-0.01; x1],[FR_prev(2); f2],'-ok')
                xlim([-0.01 0.02])
                ylim([0 2.5])
                pos_slope(av_unit)=p_all(1);
                
                
                subplot(2,4,5)
                plot(p2(1),p1(1),'ob')
                hold on
                
                
                
                subplot(2,4,6)
                plot(p2(2),p1(2),'ob')
                hold on
                
                
                
                FR1_all=simspikes(deltak(touch_idx==1)',p2,1); %first touch with later touches tuning cuve
                FR2=simspikes(deltak(touch_idx>1)',p2,1); %later touches with corresponding tuning curve
                
                xk(av_unit,:)=[mean(FR(touch_idx==1)),mean(FR1_all(:))];
                xr(av_unit,:)=[mean(FR(touch_idx==1)),mean(lambda1(:))];
                
                subplot(2,4,4)
                hold on
                plot(mean(FR(touch_idx==1)),mean(lambda1(:)),'or')
                plot(mean(FR(touch_idx==1)),mean(FR1_all(:)),'ok')
                
                
                % response to first touch is X% higher than predicted
                ratioFRb(av_unit,:)=nanmean(FR(touch_idx==1)'./FR1_all(:)');
                
                %% test if FR from first touch is higher than prediction
                [hFR(av_unit),pvalFR(av_unit)] = ttest(FR(touch_idx==1)',FR1_all(:)','Tail','right');
                %FR_decreased(av_unit)=100-100*mean(FR(touch_idx>1))/mean(FR(touch_idx==1));
                
    
                
                error1(av_unit)=mean(abs(FR(touch_idx==1)-lambda1(:)));
                error2(av_unit)=mean(abs(FR(touch_idx==1)-FR1_all(:)));

                
                
                %subplot(3,3,9)
                %hold on
                %errorbar([1 2],[mean(FR((touch_idx>1) & (whisker==1))),mean(FR((touch_idx>1) & (whisker==2)))],[std(FR((touch_idx>1) & (whisker==1))),std(FR((touch_idx>1) & (whisker==2)))])
                hw(av_unit)=ttest2(FR((whisker==1)),FR((whisker==2)));
                
                %mean(FR(touch_idx==2))
                %(mean(FR(touch_idx==1))-mean(FR(touch_idx==2)))
                %change_explained_by_adap(av_unit)=(1-ratioFR(av_unit,2))/(1-(mean(FR(touch_idx>1))/mean(FR(touch_idx==1))))
                [h_atte(av_unit),pval_atte(av_unit)]=ttest2(FR(touch_idx==1),FR(touch_idx>1));
                
                
                if h_atte(av_unit) %(mean(FR(touch_idx==1))-mean(FR(touch_idx>1)))>0
                    %change_explained_by_adap(av_unit)=min([(median(FR(touch_idx==1)-FR1_all(:)))/(mean(FR(touch_idx==1))-mean(FR(touch_idx>1))) 1e6]);
                    change_explained_by_adap(av_unit)=min([(median(FR(touch_idx==1)-FR1_all(:)))/(median(FR(touch_idx==1))-median(FR(touch_idx>1))) 1e6]);
                    
                    %change_explained_by_adap(av_unit)=(1-ratioFR(av_unit,2))/(1-(mean(FR(touch_idx>3))/mean(FR(touch_idx==1))));
                    
                else
                    change_explained_by_adap(av_unit)=NaN;
                end
                FR_touch(av_unit,:)=[mean(FR(touch_idx==1)) mean(FR(touch_idx==2)) mean(FR(touch_idx==3)) mean(FR(touch_idx>3))]/mean(FR(touch_idx==1));
                

                %
                
                
                
                %                 subplot(3,3,7)
                %                 f1_norm = polyval(p1_norm,x1_norm);
                %                 f2_norm = polyval(p2_norm,x2_norm);
                %                 plot(x1_norm,f1_norm./mean(FR(touch_idx==1)),'-or')
                %                 hold on
                %                 plot(x2_norm,f2_norm./mean(FR(touch_idx>1)),'-ok')
                %
                %                 subplot(3,3,8)
                %                 FR1_1=simspikes(deltak_norm(touch_idx==1)',p1_norm,50); %first touch with correspondent tuning curve
                %                 FR1_all=simspikes(deltak_norm(touch_idx==1)',p2_norm,50); %first touch with mixed tunning cuve
                %                 FR2=simspikes(deltak_norm(touch_idx>1)',p2_norm,50); %later touches with corresponding tuning curve
                %                 FR2_all=simspikes(deltak_norm(touch_idx>1)',p1_norm,50); %later touches with mixed tuning curve
                %                 plot(mean(FR(touch_idx==1)),mean(FR1_1(:)),'or')
                %                 hold on
                %                 plot(mean(FR(touch_idx==1)),mean(FR1_all(:)),'xr')
                %                 plot(mean(FR(touch_idx>1)),mean(FR2(:)),'ok')
                %                 plot(mean(FR(touch_idx>1)),mean(FR2_all(:)),'xk')
                %
                
                if p2(1)<p1(1)
                    counter=counter+1;
                end
                
            end
            
            clear FR FR_prev FR1_all FR2_all FR2
            
        end
    end

    
    if f==4
        figure(fig4)
        subplot(2,4,1)
        plot(sort(deltak),1:numel(deltak),'b')
        box off
        xlabel('|\Delta \kappa_{3D}| [1/mm]')
        ylabel('Touch number')
    end
    
    clear deltak touch_idx
    cd ..
end


figure(fig4)
subplot(2,4,5)
hold on
plot([-25 160],[-25 160],'-b','Linewidth',2)
plot([0 0],[-25 160],'-k','Linewidth',1)
plot([-25 160],[0 0],'-k','Linewidth',1)
box off
xlabel('Later touches slope')
ylabel('First touches slope')
ylim([-25 160])
xlim([-25 160])

plot(mean(coeff2(:,1)),mean(coeff1(:,1)),'.k')
errorbar(mean(coeff2(:,1)),mean(coeff1(:,1)),std(coeff2(:,1)),'horizontal','.b')
errorbar(mean(coeff2(:,1)),mean(coeff1(:,1)),std(coeff1(:,1)),'vertical','.b')

subplot(2,4,6)
box off
xlabel('Later touches intercept')
ylabel('First touches intercept')
ylim([0 2])
xlim([0 2])
plot([0 2],[0 2],'-b','Linewidth',2)
plot(mean(coeff2(:,2)),mean(coeff1(:,2)),'.k')
errorbar(mean(coeff2(:,2)),mean(coeff1(:,2)),std(coeff2(:,2)),'horizontal','.b')
errorbar(mean(coeff2(:,2)),mean(coeff1(:,2)),std(coeff1(:,2)),'vertical','.b')

%Does the intercept significantly decrease across the population
[hintercept,pvalintercept] = ttest(coeff1(:,2)',coeff2(:,2)');

% how much decrease
intercept_decrease=1-mean(coeff2(:,2)./coeff1(:,2));
%Does the slope significantly decrease across the population?
[hslope,pvalslope] = ttest(coeff1(:,1),coeff2(:,1));
slope_decrease=coeff2(:,1)./coeff1(:,1);

figure(fig4)
subplot(2,4,2)
xlabel('|\Delta \kappa_{3D}| [1/mm]')
ylabel('Firing rate [spikes/touch]')

subplot(2,4,3)
xlabel('|\Delta \kappa_{3D}| [1/mm]')
ylabel('Firing rate [spikes/touch]')
xticks([-0.01 0 0.01 0.02])
xticklabels({'Before touch','0','0.01','0.02'})
box off


subplot(2,4,4)
box off
xlabel('Observed FR [spikes/touch]')
ylabel('Predicted FR [spikes/touch]')
plot([0 3],[0 3],'b')
xlim([0 2.5])
ylim([0 2.5])
plot(mean(xk(:,1)),mean(xk(:,2)),'.k')
plot(mean(xr(:,1)),mean(xr(:,2)),'.r')
errorbar(mean(xk(:,1)),mean(xk(:,2)),std(xk(:,1)),'horizontal','.k')
errorbar(mean(xk(:,1)),mean(xk(:,2)),std(xk(:,2)),'vertical','.k')
errorbar(mean(xr(:,1)),mean(xr(:,2)),std(xr(:,1)),'horizontal','.r')
errorbar(mean(xr(:,1)),mean(xr(:,2)),std(xr(:,2)),'vertical','.r')


% subplot(3,3,5)
% box off
% xlim([0.8 2.2])
% ylim([0 2])
% plot([1 2],[1 1],'k')
% xticks([1 2])
% ylabel('Ratio FR tuning curves')
% xticklabels({'later touch ratio','later touches ratio'})
% errorbar([1 2], nanmean(ratioFR,1),nanstd(ratioFR,[],1),'.-k','LineWidth',2)


%nanmean(ratioFR)


% xlim([0 3])
% plot([1 2],[1 1],'k')
% ylabel('Ratio FR tuning curves')
% xticklabels({'','first touch ratio','later touches ratio',''})
%errorbar([1 2], nanmean(ratioFR_norm),nanstd(ratioFR_norm),'.-k','LineWidth',2)




% figure(fig1)
% subplot(3,3,6)
% xlabel('Touch number')
% xticks([1 2 3 4])
% ylabel('K')
% xlim([0.5 4.5])
% 
% errorbar(1:4,mean(deltak_touch,1),std(deltak_touch,[],1),'k')
% 
% [hd1,pd1] = ttest2(deltak_touch(:,1),deltak_touch(:,2));
% [hd2,pd2] = ttest2(deltak_touch(:,2),deltak_touch(:,3));
% [hd3,pd3] = ttest2(deltak_touch(:,3),deltak_touch(:,4));



% subplot(3,3,9)
% hold on
% plot([1 2 3 4],FR_touch')
% xlabel('Touch number')
% xticks([1 2 3 4])
% ylabel('FR')
% xlim([0.5 4.5])

%Are the reponses to touch significantly weaker over time?

[h1,p1] = ttest(FR_touch(:,2),1,'Tail','left');
[h2,p2] = ttest2(FR_touch(:,2),FR_touch(:,3),'Tail','right');
[h3,p3] = ttest2(FR_touch(:,3),FR_touch(:,4),'Tail','right');

% how many units show adaptation?
lower=sum(hFR)
mean(pvalFR(hFR>0))


% response to first touch is X% higher than predicted
mean(ratioFRb(ratioFRb<10))
mean(pval_atte(h_atte>0))
clear deltak touch_idx

%  figure
% histogram(change_explained_by_adap,0:0.05:1.2)
% nanmean(change_explained_by_adap)
% nanstd(change_explained_by_adap)
% sum(isinf(change_explained_by_adap))
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