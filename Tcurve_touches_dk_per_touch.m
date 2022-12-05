function [x1,x2,x3,x4,p,deltak,touch_idx,FR,whisker]=Tcurve_touches_dk_per_touch(touches,touches_whisker,psth,k_c1,k_c2,N_first)
%this code assumes that the input is dk1 and dk2
counter=1;
k_hor(:,:,1)=k_c1;
k_hor(:,:,2)=k_c2;
window=100;
window_ad=600;
%N_first=5;
N_latter=N_first;
percentile_1=1;
percentile_2=0.95;
upto=30;    
%compute prior std as the mean of dk1 and dk2
%prior_std=[std(k_c1(:)) std(k_c2(:))];

%prior_std=[std(k_hor(:)) std(k_hor(:))];

%compute prior std as std from some of the first touches (Rasmus idea)
% n_pop=10;
% select_t = floor(rand(n_pop).*size(k_c1,1))+1;
% pop_dk=[];
% for i=1:n_pop
%     idx=find(touches(select_t(i),:),1,'first');
%     w=find(diff(squeeze(touches_whisker(select_t(i),idx-1:idx,:)))==1);
%     
%     if numel(w)>1
%         continue
%     else
%         whisker=w;
%     end
%     if idx+30<=3488
%         if idx-window_ad<1
%             pop_dk=[pop_dk zeros(1,abs(idx-window_ad)) k_hor(select_t(i),1:idx+5,whisker)]; %%curvature hor
%         else
%             pop_dk=[pop_dk k_hor(select_t(i),idx-window_ad:idx+5,whisker)]; %%curvature hor
%         end
%     end
% end
% 
% prior_std=[std(abs(pop_dk)) std(abs(pop_dk))];
% 
% clear whisker idx
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
            
            if idx(t)-window_ad<1
                        std_dk(counter)=std(abs(k_hor(trial,idx(t)-500:idx(t)-1,whisker(counter)))); %%curvature hor
                        mean_dk(counter)=std(abs(k_hor(trial,idx(t)-500:idx(t)-1,whisker(counter))));
            else
                        std_dk(counter)=std(abs(k_hor(trial,idx(t)-window_ad:idx(t)-1,whisker(counter))));
                        mean_dk(counter)=mean(abs(k_hor(trial,idx(t)-window_ad:idx(t)-1,whisker(counter))));
            end
                    %std_dk(counter)=std(sum(k_hor(trial,idx(t)-window_ad:idx(t)-1,:),3));
                    %std_dk(counter)=std(k_hor(trial,idx(t)-window_ad:idx(t)-1,whisker(counter)));
                
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
% std_dk=std_dk(whisker==1);

FR=sum(spikes(:,100:100+upto),2);

deltak=abs(deltak_hor);
%first touch
FR_1=sum(spikes(touch_idx==1,100:100+upto),2);
deltak_1=abs(deltak_hor(touch_idx==1));
%for the first touch, std of delta kappa should be zero by definition

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
    p1=nan;
    p2=nan;
    return
end
  
    [~,x1,p(1,:)]=plotTCurve(deltak_1,FR_1,'',N_first,percentile_1);
    [~,x2,p(2,:)]=plotTCurve(deltak_2,FR_2,'',N_latter,percentile_2);
    [~,x3,p(3,:)]=plotTCurve(deltak_3,FR_3,'',N_latter,percentile_2);
    [~,x4,p(4,:)]=plotTCurve(deltak_4,FR_4,'',N_latter,percentile_2);
end

function [Nsample,xtemp,P,y]=plotTCurve(Amp,Nspikes,plot_d,Nbins,percentile)
Nspikes(isinf(Amp))=[];
Amp(isinf(Amp))=[];
% [~,idx_temp]=sort(Amp);
% plot(Nspikes(idx_temp))
if Nbins<2
    Nsample=nan;
    xtemp=nan;
    P=nan;
    return
end
%linspace(0,0.95,Nbins+1)
x = quantile(Amp,linspace(0,percentile,Nbins+1));
%max(x)
% Amin=min(Amp);
% Amax=max(Amp);
% x=linspace(Amin,Amax,Nbins+1);
%x=[0.5 0.8 0.9 0.95 1 1.05 1.1 1.2 1.4 2];
y2=discretize(Amp,x)';
y=zeros(Nbins,1);
err=y;
Nsample=y;
ind=1;
for j=1:Nbins
    
    idx=find(y2==j);
    if ~isempty(idx)
        %Firing probability
        %y(ind)=sum(Nspikes(idx)>0)/numel(Nspikes(idx));
        %FR

        y(ind)=mean(Nspikes(idx));
        err(ind)=std(Nspikes(idx));
        Nsample(ind)=numel(Nspikes(idx));
    end
    ind=ind+1;
end
%y=y/(sum(Nspikes)/numel(Nspikes));
xtemp=(x(1:end-1)+x(2:end))'/2;
P=polyfit(xtemp,y,1);
f = polyval(P,xtemp);


if strcmp(plot_d,'')
else
    
    errorbar(xtemp,y,err,[ '.' plot_d])
    hold on
    plot(xtemp,f,['-o' plot_d])
    title(num2str(P(1)))
end
end