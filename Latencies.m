function [latencies,touch_idx,FR_prev_tmp]=Latencies(touches,touches_whisker,psth,do_plot,subplot_idx)
%this code assumes that the input is dk1 and dk2
counter=1;
upto=50;

latencies=[];
touch_idx=[];

FR_prev_tmp=[];

Ntrials=size(touches,1);

for trial=1:Ntrials
    idx=find(touches(trial,:));
    for t=1:size(idx,2)
        if t>=2 && ((idx(t)-idx(1))<=20)
            continue
        else
            
            
            if idx(t)+upto<=3488
                first_sp=find(psth(trial,idx(t)+1:idx(t)+upto),1,'First');

                if ~isempty(first_sp)
                
                    latencies(counter)=first_sp;

                    touch_idx(counter)=t;%% order by trial

                    FR_prev_tmp(counter)=sum(psth(trial,idx(t)-50:idx(t)-1))*1000/50;
                end
               
            end
                  
            counter=counter+1;
        end
    end
    
end

%% first touch
latencies_touch=[mean(latencies(touch_idx==1)) mean(latencies(touch_idx==2)) mean(latencies(touch_idx==3)) mean(latencies(touch_idx>=4))];
latencies_touch_std=[std(latencies(touch_idx==1)) std(latencies(touch_idx==2)) std(latencies(touch_idx==3)) std(latencies(touch_idx>=4))];


if do_plot
       
    subplot(4,6,subplot_idx)
    hold on
    errorbar(1:4,latencies_touch(1:4),latencies_touch_std(1:4),'Color',[0.5 0.5 0.5])
    xlabel('Touch Number')
    ylabel('Latencies [ms]')
    box off 

    
end

end