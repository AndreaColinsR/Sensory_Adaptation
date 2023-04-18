function [x1,x2,x3,x4,p,deltak,touch_idx,FR,whisker]=Tcurve_touches_dk_per_touch(touches,touches_whisker,psth,k_c1,k_c2,NpointsTcurve)
%% Tcurve_touches_dk_per_touch computes 4 tuning curves: one for only first, only second and only third touches of the sequences and one for touches >4
% INPUTS
% 
% touches: binary matrix, rows are trials, columns are time bins in ms. 
%     1= at least 1 whisker is touching. 
%     0= no whiskers are touching. 
% 
% touches_whisker:binary 3D matrix [trials, timebins,whiskers] 1 indicates
% that a whisker is touching the pole.
% 
% psth: binary matrix, rows are trials, columns are time bins in ms. 1
% indicate a spike.
% 
% k_c1= Delta kappa of whisker 1. [trials, timebins]
%
% k_c2= Delta kappa of whisker 2. [trials, timebins]
% 
% NpointsTcurve= number of bins to segment the data for the tuning curve
%
% OUTPUTS
%
% x1,2,3,4= medians of the bins of delta k for the four tuning curves
%
% p= [tuning curve number,2]= the first column is the slope and the second
% is the intercept of the tuning curve
% 
% deltak= delta kappa values for all touches in the session
%
% touch_idx= number of each touch in the trial
%
% FR= Number of spikes from the touch onset upto 30 ms after
%
% whisker= whisker identity of the touch.



%this code assumes that the input is dk1 and dk2
counter=1;
k_hor(:,:,1)=k_c1;
k_hor(:,:,2)=k_c2;
window=100;
percentile_1=1;
percentile_2=0.95;
upto=30;    

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
% std_dk=std_dk(whisker==1);

FR=sum(spikes(:,100:100+upto),2);

deltak=abs(deltak_hor);
%first touch
FR_1=sum(spikes(touch_idx==1,100:100+upto),2);
deltak_1=abs(deltak_hor(touch_idx==1));

%latter touches
FR_2=sum(spikes(touch_idx==2,100:100+upto),2);
deltak_2=abs(deltak_hor(touch_idx==2));

FR_3=sum(spikes(touch_idx==3,100:100+upto),2);
deltak_3=abs(deltak_hor(touch_idx==3));

FR_4=sum(spikes(touch_idx>4,100:100+upto),2);
deltak_4=abs(deltak_hor(touch_idx>4));

FR2=sum(spikes(touch_idx==2,100:130),2);
total_sp=sum(FR_1)+sum(FR2);
if total_sp<20
    x1=nan;
    x2=nan;
    return
end
  
    [~,x1,p(1,:)]=plotTCurve(deltak_1,FR_1,'',NpointsTcurve,percentile_1);
    [~,x2,p(2,:)]=plotTCurve(deltak_2,FR_2,'',NpointsTcurve,percentile_2);
    [~,x3,p(3,:)]=plotTCurve(deltak_3,FR_3,'',NpointsTcurve,percentile_2);
    [~,x4,p(4,:)]=plotTCurve(deltak_4,FR_4,'',NpointsTcurve,percentile_2);
end