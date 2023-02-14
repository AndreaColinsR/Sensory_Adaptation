function [hit,hit3,fa,fa2,fa3]=Detection_test(Data,total_psth,varargin)
%aim: compare detectability of latter touches with dectectability of

idx=ones(size(total_psth,1),1)>0;
 
if nargin==3
    idx=varargin{1};
end
   
% evaluate histogram of spike counts
k=0:1:21;
nk=numel(k)-1;
hit=nan(1,nk);
fa=nan(1,nk);
hit3=nan(1,nk);
fa2=nan(1,nk);
fa3=nan(1,nk);


% load('all_period_touches.mat','touches_matrix')
% period_touches=touches_matrix(idx,:,:);
% load('all_touches.mat','touches_matrix')
% touches_matrix=touches_matrix(idx,:,:);

period_touches=Data.all_period_touch(idx,:,:);

touches_matrix=Data.touch(idx,:,:);


counter=1;
counter2=1;

for i=1:size(total_psth,1)
    idx=find(touches_matrix(i,:));
    
    %%%% first touch
    if numel(idx)>=2
        if (idx(2)-idx(1))>20
            after_touch(counter,:)=total_psth(i,idx(1):idx(1)+30);
            counter=1+counter;
        end
        
        %%% later touches
        for t=2:size(idx,2)

            if (idx(t)-idx(1))>20
                if idx(t)+30<=3488
                    other_touch(counter2,:)=total_psth(i,idx(t):idx(t)+30);
                    between_t=total_psth(i,idx(t)-35:idx(t)-5);
                    FR4(counter2)=sum(between_t);
                    counter2=counter2+1;
                end
            end
        end
        
    end
end


FR2=sum(after_touch,2);
FR3=sum(other_touch,2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% before first touch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noise=total_psth(:,1:500);

counter=1;
window=30;
bins=500/window;

for j=1:size(noise,1)
    for i=1:bins
        FR(counter)=sum(noise(j,window*(i-1)+1:window*(i-1)+window),2);
        counter=1+counter;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counter=1;

clear noise
for j=1:size(total_psth,1)
    idx=find(touches_matrix(j,:),1,'first');
    
    noise=total_psth(j,idx+30:1900);
    tmp_touch=period_touches(j,idx+30:1900);
    if ~isempty(noise)
        bins=floor(size(noise,2)/window);
        for i=1:bins
            if sum(tmp_touch(window*(i-1)+1:window*(i-1)+window))==0
                FR6(counter)=sum(noise(window*(i-1)+1:window*(i-1)+window),2);
                counter=1+counter;
            end
        end
        clear noise
    end
end

pois_1=histcounts(FR,-0.5:1:21.5,'normalization','probability');
pois_2=histcounts(FR2,-0.5:1:21.5,'normalization','probability');
pois_3=histcounts(FR3,-0.5:1:21.5,'normalization','probability');
pois_4=histcounts(FR4,-0.5:1:21.5,'normalization','probability');
pois_6=histcounts(FR6,-0.5:1:21.5,'normalization','probability');



for i=1:size(k,2)-1
    hit(i)=trapz(k(i:end),pois_2(i:end));
    hit3(i)=trapz(k(i:end),pois_3(i:end));
    fa(i)=trapz(k(i:end),pois_1(i:end));
    fa2(i)=trapz(k(i:end),pois_4(i:end));
    fa3(i)=trapz(k(i:end),pois_6(i:end));
end


end