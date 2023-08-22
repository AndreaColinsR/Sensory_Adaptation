function [change_explained_by_adap,FR_touch]=all_Tcurves(fig4)
do_extra_plot=0;
counter=0;
counter2=1;
av_unit=0;
percentile=0.95;
NpointsTcurve=4;
ff=dir('*.mat*');

nunits=33;
FR_touch=nan(nunits,4);
h_atte=nan(nunits,1);
coeff_first=nan(nunits,2);
coeff_later=nan(nunits,2);
hw=nan(nunits,1);
pval_atte=nan(nunits,1);
pvalFR=nan(nunits,1);
hFR=nan(nunits,1);
xk=nan(nunits,2);
xr=nan(nunits,2);
change_explained_by_adap=nan(nunits,1);

%% Define population FR per mouse
mouse_id=zeros(size(ff,1),1);

for f=1:size(ff,1)
    mouse_id(f)=str2double(ff(f).name(4:5));
end

mouse_number=unique(mouse_id);
intercept_mouse=cell(numel(mouse_number),1);

for f=1:size(ff,1)
    
    load(ff(f).name,'Data')
    this_mouse=str2double(ff(f).name(4:5))==mouse_number;
    
    
    %total_psth=zeros(size(Data.touch));
    for i_unit=1:size(Data.unit,2)
        
        total_psth=Data.unit(i_unit).spikes;
        
        
        delete=size(total_psth,2)-size(Data.deltak_w1,2);
        Data.deltak_w1=[zeros(size(Data.deltak_w1,1),delete),Data.deltak_w1];
        Data.deltak_w2=[zeros(size(Data.deltak_w1,1),delete),Data.deltak_w2];
        
        clear deltak touch_idx FR
        
        [x1,p1,~,p2,~,~,deltak,touch_idx,FR,~,FR_prev,whisker]=Tcurve_touches_dk(Data.touch,Data.touch_per_whisker,total_psth,Data.deltak_w1,Data.deltak_w2,do_extra_plot,NpointsTcurve,0,percentile,percentile);
        deltak=abs(deltak);
        
        [lambda1,~]=Tcurve_touches_dk_cross_val(Data.touch,Data.touch_per_whisker,total_psth,Data.deltak_w1,Data.deltak_w2,0,NpointsTcurve);
        
        if ~isnan(x1)
            av_unit=av_unit+1;
            
            %% testing whisker preference
            % does this neuron show a preference for a whisker
            hw(av_unit)=ttest2(FR((whisker==1 & touch_idx>4)),FR((whisker==2 & touch_idx>4)));
            
%             if hw(av_unit)
%                 subplot(3,3,9)
%                 hold on
%                 errorbar([1 2],[mean(FR((touch_idx==1) & (whisker==1))),mean(FR((touch_idx==1) & (whisker==2)))],[std(FR((touch_idx==1) & (whisker==1))),std(FR((touch_idx==1) & (whisker==2)))])
%                 
%                 %choose the prefered whisker
%                 %[~,pref_w]=max([mean(FR((whisker==1 & touch_idx>4))),mean(FR((whisker==2 & touch_idx>4)))])
%                 
%                 % choose the whisker that touches more times
%                 [~,pref_w]=max([sum(whisker==1 & touch_idx>4),sum((whisker==2 & touch_idx>4))]);
%                 
%                 
%                 selected_touch=touch_idx==1 & whisker==pref_w;
%                 nsamples_1=sum(selected_touch);
%                 [~,x1,p1]=plotTCurve(deltak(selected_touch),FR(selected_touch),'',NpointsTcurve,percentile);
%                 
%                 selected_touch=touch_idx>1 & whisker==pref_w;
%                 nsamples_2=sum(selected_touch);
%                 [~,~,p2]=plotTCurve(deltak(selected_touch),FR(selected_touch),'',NpointsTcurve,percentile);
%                 
%                 FR=FR(whisker==pref_w);
%                 deltak=deltak(whisker==pref_w);
%                 touch_idx=touch_idx(whisker==pref_w);
%                 %whisker=whisker(whisker==pref_w);
%                 
%             end
            
            coeff_first(av_unit,:)=p1;
            coeff_later(av_unit,:)=p2;
            
            f1 = polyval(p1,x1);
            f2 = polyval(p2,x1);
            
            figure(fig4)
            if strcmp('Glu35_10112017H.mat',ff(f).name) && counter2==1
                
                subplot(2,4,2)
                nsamples=plotTCurve(deltak(touch_idx>1),FR(touch_idx>1),'k',NpointsTcurve,percentile);
                nsamples=plotTCurve(deltak(touch_idx==1),FR(touch_idx==1),'r',NpointsTcurve,percentile);
                
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
            
            % slope
            subplot(2,4,5)
            plot(p2(1),p1(1),'ob')
            hold on
            
            
            % intercept
            subplot(2,4,6)
            plot(p2(2),p1(2),'ob')
            hold on
            
            intercept_mouse{this_mouse}=[intercept_mouse{this_mouse}; [p1(2) p2(2)]];
            
            
            % check here for preference
            FR1_all=simspikes(deltak(touch_idx==1)',p2,1); %first touch with later touches tuning cuve
            % check here for preference
            xk(av_unit,:)=[mean(FR(touch_idx==1)),mean(FR1_all(:))];
            xr(av_unit,:)=[mean(FR(touch_idx==1)),mean(lambda1(:))];
            
            % check here for preference
            subplot(2,4,4)
            hold on
            plot(mean(FR(touch_idx==1)),mean(lambda1(:)),'or')
            plot(mean(FR(touch_idx==1)),mean(FR1_all(:)),'ok')
            
            
            % response to first touch is X% higher than predicted
            %ratioFRb(av_unit)=mean(FR(touch_idx==1)'./FR1_all(:)','omitnan');
            
            %% test if FR from first touch is higher than prediction
            [hFR(av_unit),pvalFR(av_unit)] = ttest(FR(touch_idx==1)',FR1_all(:)','Tail','right');
      
            
            
            %% t-test first vs later responses
            
            [h_atte(av_unit),pval_atte(av_unit)]=ttest2(FR(touch_idx==1),FR(touch_idx>1));
            
            
            if h_atte(av_unit) %(mean(FR(touch_idx==1))-mean(FR(touch_idx>1)))>0
                %change_explained_by_adap(av_unit)=min([(median(FR(touch_idx==1)-FR1_all(:)))/(mean(FR(touch_idx==1))-mean(FR(touch_idx>1))) 1e6]);
                change_explained_by_adap(av_unit)=min([(median(FR(touch_idx==1)-FR1_all(:)))/(median(FR(touch_idx==1))-median(FR(touch_idx>1))) 1e6]);
                
                %change_explained_by_adap(av_unit)=(1-ratioFR(av_unit,2))/(1-(mean(FR(touch_idx>3))/mean(FR(touch_idx==1))));
                
            else
                change_explained_by_adap(av_unit)=NaN;
            end
            FR_touch(av_unit,:)=[mean(FR(touch_idx==1)) mean(FR(touch_idx==2)) mean(FR(touch_idx==3)) mean(FR(touch_idx>3))]/mean(FR(touch_idx==1));
            
            
            
            
            if p2(1)<p1(1)
                counter=counter+1;
            end
            
        end
        
        clear FR FR_prev FR1_all FR2_all FR2
        
    end
    
    
    if strcmp('Glu35_10112017H.mat',ff(f).name)
        figure(fig4)
        subplot(2,4,1)
        plot(sort(deltak),1:numel(deltak),'b')
        box off
        xlabel('|\Delta \kappa_{3D}| [1/mm]')
        ylabel('Touch number')
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

plot(mean(coeff_later(:,1)),mean(coeff_first(:,1)),'.k')
errorbar(mean(coeff_later(:,1)),mean(coeff_first(:,1)),std(coeff_later(:,1)),'horizontal','.b')
errorbar(mean(coeff_later(:,1)),mean(coeff_first(:,1)),std(coeff_first(:,1)),'vertical','.b')

subplot(2,4,6)
box off
xlabel('Later touches intercept')
ylabel('First touches intercept')
ylim([0 2])
xlim([0 2])
plot([0 2],[0 2],'-b','Linewidth',2)
plot(mean(coeff_later(:,2)),mean(coeff_first(:,2)),'.k')
errorbar(mean(coeff_later(:,2)),mean(coeff_first(:,2)),std(coeff_later(:,2)),'horizontal','.b')
errorbar(mean(coeff_later(:,2)),mean(coeff_first(:,2)),std(coeff_first(:,2)),'vertical','.b')

%Does the intercept significantly decrease across the population
%[hintercept,pvalintercept] = ttest(coeff1(:,2)',coeff2(:,2)');

% how much decrease
% intercept_decrease=1-mean(coeff2(:,2)./coeff1(:,2));
% Does the slope significantly decrease across the population?
% [hslope,pvalslope] = ttest(coeff1(:,1),coeff2(:,1));
% slope_decrease=coeff2(:,1)./coeff1(:,1);

disp('----------------------------------------------')
disp(['Proportion of units showing positive slope =' num2str(sum(coeff_later(:,1)>0)/nunits)])
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

disp('How many units do show adaptation?')
disp(num2str(sum(hFR,'omitnan')))
disp('-------------------------------')


disp('')
disp('')
disp('-----------------------------------------')
disp('Intercept decrease per mouse (later/first)')
disp(['Whole population  = ' num2str(mean(coeff_later(:,2)./coeff_first(:,2),'omitnan'))])

for im=1:numel(mouse_number)
    
    disp(['Mouse ' num2str(im) ' =' num2str(mean(intercept_mouse{im}(:,2)./intercept_mouse{im}(:,1),1))])
    
end

%% Whisker preference
sum(hw,'omitnan')
end

function FR=simspikes(x,p,Nsim)
lambda=p(1)*x+p(2);
lambda=repmat(lambda,1,Nsim);
FR = lambda;
end