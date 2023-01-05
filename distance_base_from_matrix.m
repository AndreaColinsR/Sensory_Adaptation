function [distance_1,distance_2,distance_3,distance_4,distance,touch_idx,hit]=distance_base_from_matrix(idx_correct)
touch_idx=[];
hit=[];
counter=1;
info=xlsread('ledtrials.xlsx');
go_trials=find(info(:,7)==2);
idx_correct_1=info(go_trials,6)==1;

load('touches_whisker.mat','touches_whisker')
load('distance_matrix.mat','distance_matrix')


if numel(go_trials)>sum(idx_correct)
    touches_whisker=touches_whisker(idx_correct,:,:);
    distance_matrix=distance_matrix(idx_correct,:,:);
end


for trial=1:sum(idx_correct)
    
    %find touch onset
    idx1=find(diff(touches_whisker(trial,:,1))>0);
    idx2=find(diff(touches_whisker(trial,:,2))>0);
    [idx,ind]=sort([idx1 idx2]);
    whisker_trial=[ones(size(idx1)) ones(size(idx2))*2];
    whisker_trial=whisker_trial(ind);
    
    
    
    for t=1:length(whisker_trial)
        if idx(t)<3478
            
            distance(counter)=distance_matrix(trial,idx(t)+1,whisker_trial(t));
            touch_idx(counter)=t;
            hit(counter)=idx_correct_1(trial);
            counter=counter+1;
            
        end
    end
    
end
distance_1=distance(touch_idx==1);
distance_2=distance(touch_idx==2);
distance_3=distance(touch_idx==3);
distance_4=distance(touch_idx>3);

end