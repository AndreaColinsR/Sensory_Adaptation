function [change_explained_by_adap,FR_touch]=all_Tcurves_per_whisker(fig4)
do_extra_plot=0;
counter=0;
counter2=1;
av_unit=0;
percentile_1=0.95;
percentile_2=0.95;
NpointsTcurve=4;
ff=dir('*.mat*');

nunits=33;
FR_touch=nan(nunits,4);
h_atte=nan(nunits,1);
coeff1=nan(nunits,2);
coeff2=nan(nunits,2);

for f=1:size(ff,1)
    
    load(ff(f).name,'Data')
    
    total_psth=zeros(size(Data.touch));
    for i_unit=1:size(Data.unit,2)
        
        total_psth=Data.unit(i_unit).spikes;
        
        
        delete=size(total_psth,2)-size(Data.deltak_w1,2);
        Data.deltak_w1=[zeros(size(Data.deltak_w1,1),delete),Data.deltak_w1];
        Data.deltak_w2=[zeros(size(Data.deltak_w1,1),delete),Data.deltak_w2];
        
        clear deltak touch_idx FR
        
        [x1,p1,~,p2,x1_norm,p1_norm,x2_norm,p2_norm,x_all,p_all,deltak,touch_idx,FR,~,FR_prev,whisker]=Tcurve_touches_dk(Data.touch,Data.touch_per_whisker,total_psth,Data.deltak_w1,Data.deltak_w2,do_extra_plot,NpointsTcurve,0,percentile_1,percentile_2);
        deltak=abs(deltak);
        
        %subplot(6,6,av_unit+1)
        [~,~,P1,~]=plotTCurve(deltak(whisker==1 & touch_idx>1),FR(whisker==1 & touch_idx>1),'',NpointsTcurve,percentile_1);
        [~,~,P2,~]=plotTCurve(deltak(whisker==2 & touch_idx>1),FR(whisker==2 & touch_idx>1),'',NpointsTcurve,percentile_1);
        
        if ~isnan(x1)
            av_unit=av_unit+1;
          [h,p]=ttest2(FR(whisker==1 & touch_idx>1),FR(whisker==2 & touch_idx>1));

           plot(av_unit,h,'o')
           hold on
        end
        
        clear FR FR_prev FR1_all FR2_all FR2
        
    end
  
    clear deltak touch_idx
    
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
%[hintercept,pvalintercept] = ttest(coeff1(:,2)',coeff2(:,2)');

% how much decrease
% intercept_decrease=1-mean(coeff2(:,2)./coeff1(:,2));
% Does the slope significantly decrease across the population?
% [hslope,pvalslope] = ttest(coeff1(:,1),coeff2(:,1));
% slope_decrease=coeff2(:,1)./coeff1(:,1);

disp('----------------------------------------------')
disp(['Proportion of units showing positive slope =' num2str(sum(coeff2(:,1)>0)/nunits)])
disp('----------------------------------------------')
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


%Are the reponses to touch significantly weaker over time?
% 
% [h1,p1] = ttest(FR_touch(:,2),1,'Tail','left');
% [h2,p2] = ttest2(FR_touch(:,2),FR_touch(:,3),'Tail','right');
% [h3,p3] = ttest2(FR_touch(:,3),FR_touch(:,4),'Tail','right');

% how many units show adaptation?
lower=sum(hFR)
mean(pvalFR(hFR>0))


% response to first touch is X% higher than predicted
% mean(ratioFRb(ratioFRb<10))
% mean(pval_atte(h_atte>0))
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