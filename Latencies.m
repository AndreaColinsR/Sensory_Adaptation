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
            
            w=find(diff(squeeze(touches_whisker(trial,idx(t)-1:idx(t),:)))==1);
            
            if numel(w)>1
                continue
            else
                whisker(counter)=w;
            end
            
            
            if idx(t)+upto<=3488
                first_sp=find(psth(trial,idx(t)+1:idx(t)+upto),1,'First');

                if ~isempty(first_sp)
                
                    latencies(counter)=first_sp;

                    touch_idx(counter)=t;%% order by trial

                    FR_prev_tmp(counter)=sum(psth(trial,idx(t)-50:idx(t)-1))*1000/50;
                    %% to do: peak latency. 
                end
                
                

            end
            
            
            
            
            
            % [end2,a]=min([idx(t)+window size(psth,2)]);
            % if a==1
            %     spikes(counter,:)=psth(trial,idx(t)-100:idx(t)+window);
            % else
            %     aux=zeros(1,size(spikes(1,:),2)-size(psth(trial,idx(t)-100:end2),2));
            %     spikes(counter,:)=[psth(trial,idx(t)-100:end2) aux];
            % 
            % end
            counter=counter+1;
        end
    end
    
end

%% first touch
%FR_prev_1=mean(FR_prev_tmp(touch_idx==1));
%FR_1=sum(spikes(touch_idx==1,100:100+upto),2);
latencies_touch=[mean(latencies(touch_idx==1)) mean(latencies(touch_idx==2)) mean(latencies(touch_idx==3)) mean(latencies(touch_idx>=4))];
latencies_touch_std=[std(latencies(touch_idx==1)) std(latencies(touch_idx==2)) std(latencies(touch_idx==3)) std(latencies(touch_idx>=4))];



%% ttest 

%latter touches
%FR_2=sum(spikes(touch_idx>1,100:100+upto),2);
% FR_prev_2=mean(FR_prev_tmp(touch_idx>1));
% deltak_2=abs(latencies(touch_idx>1));
% 
% FR_prev=[FR_prev_1 FR_prev_2];
% 
% 
% FR2=sum(spikes(touch_idx==1,100:130),2);
% total_sp=sum(FR_1)+sum(FR2);



if do_plot
       
    subplot(7,5,subplot_idx)
    hold on
    errorbar([1:4],latencies_touch([1:4]),latencies_touch_std([1:4]),'Color',[0.5 0.5 0.5])
    xlabel('Touch Number')
    ylabel('Latencies [ms]')
    box off 

    % hold on
    % plot(latencies,FR_prev_tmp,'.')
    % xlabel('Touch Number')
    % ylabel('Latencies [ms]')
    % box off 
    
    %disp([' corr latency vs FR pre touch =' num2str(corr(latencies',FR_prev_tmp'))])
    
end

end