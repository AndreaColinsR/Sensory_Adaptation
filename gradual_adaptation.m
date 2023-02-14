function gradual_adaptation(fig4)
figure(fig4)
av_unit=0;
ff=dir('*.mat*');

for f=1:size(ff,1)
    

    load(ff(f).name,'Data')
    
    %total_psth=zeros(size(Data.touch));
    for i_unit=1:size(Data.unit,2)
        
            total_psth=Data.unit(i_unit).spikes;
            
            delete=size(total_psth,2)-size(Data.deltak_w1,2);
            Data.deltak_w1=[zeros(size(Data.deltak_w1,1),delete),Data.deltak_w1];
            Data.deltak_w2=[zeros(size(Data.deltak_w1,1),delete),Data.deltak_w2];
            
            clear deltak touch_idx FR
            
            [~,~,~,~,p,deltak,touch_idx,FR]=Tcurve_touches_dk_per_touch(Data.touch,Data.touch_per_whisker,total_psth,Data.deltak_w1,Data.deltak_w2,4);
            [h_atte,~]=ttest2(FR(touch_idx==1),FR(touch_idx>1));
            
            if h_atte
                av_unit=av_unit+1;

                intercept(av_unit,:)=p(:,2)';
                
                subplot(2,4,8)
                plot(p(:,2),'-','Color',[0.5 0.5 0.5])
                hold on

                
                FR1_all=simspikes(deltak(touch_idx==1)',p(4,:),1); %predict first touch responses with later curves
                FR2_all=simspikes(deltak(touch_idx==2)',p(4,:),1); %predict second touch responses with later curves
                FR3_all=simspikes(deltak(touch_idx==3)',p(4,:),1); %predict second touch responses with later curves
                FR4_all=simspikes(deltak(touch_idx==4)',p(4,:),1); %predict second touch responses with later curves
                
                mean_pred=[mean(FR1_all(:)) mean(FR2_all(:)) mean(FR3_all(:)) mean(FR4_all(:))];
                mean_FR=[mean(FR(touch_idx==1)), mean(FR(touch_idx==2)),mean(FR(touch_idx==3)),mean(FR(touch_idx==4))];
                ratioFR(av_unit,:)=mean_pred./mean_FR;
                
                subplot(2,4,7)
                plot(mean_pred./mean_FR,'-','Color',[0.5 0.5 0.5])
                hold on
                
                clear FR p deltak touch_idx whisker
            end

    end

    
end
subplot(2,4,7)
errorbar(1:4,mean(ratioFR),std(ratioFR),'.-k','LineWidth',2)
xlim([0.5 4.5])
xlabel('Touch number')
ylabel('Predicted FR/ Observed FR')
box off
xlim([0.5 4.5])

subplot(2,4,8)
errorbar(1:4,mean(intercept),std(intercept),'.-k','LineWidth',2)
xlabel('Touch number')
ylabel('Intercept')
xlim([0.5 4.5])
box off 

end

function FR=simspikes(x,p,Nsim)
lambda=p(1)*x+p(2);
lambda=repmat(lambda,1,Nsim);
FR = lambda;
end