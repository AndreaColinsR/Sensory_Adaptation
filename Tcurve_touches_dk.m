function [x1,p1,x2,p2,x_all,p_all,deltak,touch_idx,FR,deltak_norm,FR_prev,whisker,spikes]=Tcurve_touches_dk(touches,touches_whisker,psth,k_c1,k_c2,do_plot,NpointsTcurve,subplot_idx,percentile_1,percentile_2)
%this code assumes that the input is dk1 and dk2
counter=1;
k_hor(:,:,1)=k_c1;
k_hor(:,:,2)=k_c2;
window=100;
upto=30;

%counter2=1;
% if do_plot
%     figure
%
% end


for trial=1:size(k_c1,1)
    idx=find(touches(trial,:));
    for t=1:size(idx,2)
        if t>=2 && ((idx(t)-idx(1))<=20)
            continue
        else
            
            w=find(diff(squeeze(touches_whisker(trial,idx(t)-1:idx(t),:)))==1);
            
            if numel(w)>1
                continue
            else
                whisker(counter)=w;
            end
            
            
            if idx(t)+upto<=3488
                
                deltak_hor(counter)=mean(k_hor(trial,idx(t):idx(t)+5,whisker(counter)));%%curvature hor
                
            else
                deltak_hor(counter)=mean(k_hor(trial,idx(t):end,whisker(counter)));%%curvature hor
            end
            
            touch_idx(counter)=t;%% order by trial
            
            
            FR_prev_tmp(counter)=sum(psth(trial,idx(t)-500:idx(t)-1))*30/500;
            
            
            [end2,a]=min([idx(t)+window size(psth,2)]);
            if a==1
                spikes(counter,:)=psth(trial,idx(t)-100:idx(t)+window);
            else
                aux=zeros(1,size(spikes(1,:),2)-size(psth(trial,idx(t)-100:end2),2));
                spikes(counter,:)=[psth(trial,idx(t)-100:end2) aux];
                
            end
            counter=counter+1;
        end
    end
    
end
%all together
% spikes=spikes(whisker==1,:);
% deltak_hor=deltak_hor(whisker==1);
% touch_idx=touch_idx(whisker==1);
FR=sum(spikes(:,100:100+upto),2);
deltak=deltak_hor;
deltak_norm=abs(deltak_hor);


%% first touch
FR_prev_1=mean(FR_prev_tmp(touch_idx==1));
FR_1=sum(spikes(touch_idx==1,100:100+upto),2);
deltak_1=abs(deltak_hor(touch_idx==1));


%latter touches
FR_2=sum(spikes(touch_idx>1,100:100+upto),2);
FR_prev_2=mean(FR_prev_tmp(touch_idx>1));
deltak_2=abs(deltak_hor(touch_idx>1));

FR_prev=[FR_prev_1 FR_prev_2];


FR2=sum(spikes(touch_idx==1,100:130),2);
total_sp=sum(FR_1)+sum(FR2);

if total_sp<40
    x1=nan;
    x2=nan;
    p1=nan;
    p2=nan;
    x_all=nan;
    p_all=nan;
    return
end


if do_plot
    
    
    subplot(7,5,subplot_idx)
    hold on
    [~,x2,p2]=plotTCurve(deltak_2,FR_2,'k',NpointsTcurve,percentile_2);
    [~,x_all,p_all]=plotTCurve(deltak,FR,'',NpointsTcurve,percentile_2);
    [~,x1,p1]=plotTCurve(deltak_1,FR_1,'r',NpointsTcurve,percentile_1);
    xlabel('\Delta \kappa')
    ylabel('FR')
    box off
else
    [~,x1,p1]=plotTCurve(deltak_1,FR_1,'',NpointsTcurve,percentile_1);
    [~,x2,p2]=plotTCurve(deltak_2,FR_2,'',NpointsTcurve,percentile_2);
    [~,x_all,p_all]=plotTCurve(deltak,FR,'',NpointsTcurve,percentile_1);
    
    
end

end