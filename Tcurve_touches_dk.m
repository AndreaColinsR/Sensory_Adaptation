function [x1,p1,x2,p2,x1_norm,p1_norm,x2_norm,p2_norm,x_all,p_all,deltak,touch_idx,FR,deltak_norm,FR_prev,whisker]=Tcurve_touches_dk(touches,touches_whisker,psth,k_c1,k_c2,do_plot,N_first)
%this code assumes that the input is dk1 and dk2
counter=1;
k_hor(:,:,1)=k_c1;
k_hor(:,:,2)=k_c2;
window=100;
window_ad=600;
normalised=0;
%N_first=5;
N_latter=N_first;
percentile_1=1;
percentile_2=0.95;
%counter2=1;
if do_plot 
    figure

end
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
deltak_norm=abs(deltak_hor)./std_dk;
%first touch
FR_prev_1=mean(FR_prev_tmp(touch_idx==1));
FR_1=sum(spikes(touch_idx==1,100:100+upto),2);
FR_1wong=sqrt(sum(spikes(touch_idx==1,100:100+upto),2).*FR_prev_1);
deltak_1=abs(deltak_hor(touch_idx==1));
deltak_1_norm=deltak_norm(touch_idx==1);
%FR_prev_1=mean(sum(spikes(touch_idx==1,100-upto:100-1),2));
%for the first touch, std of delta kappa should be zero by definition

%latter touches
FR_2=sum(spikes(touch_idx>1,100:100+upto),2);
%FR_prev_2=mean(sum(spikes(touch_idx>1,100-upto:100-1),2));
FR_prev_2=mean(FR_prev_tmp(touch_idx>1));
deltak_2=abs(deltak_hor(touch_idx>1));
%deltak_2_norm=abs(deltak_hor(touch_idx>1))./std_dk(touch_idx>1);
deltak_2_norm=deltak_norm(touch_idx>1);

FR_prev=[FR_prev_1 FR_prev_2];
%delete outliers of deltak
% idx_out=find(deltak_2>max(deltak_1));
% FR_2(idx_out)=[];
% deltak_2(idx_out)=[];
% 
% idx_out2=find(deltak>max(deltak_1));
% FR(idx_out2)=[];
% deltak(idx_out2)=[];
% 
% deltak_2_norm(idx_out)=[];
% touch_idx(idx_out2)=[];
% deltak_hor(idx_out2)=[];

FR2=sum(spikes(touch_idx==1,100:130),2);
total_sp=sum(FR_1)+sum(FR2);
if total_sp<40
    x1=nan;
    x2=nan;
    p1=nan;
    p2=nan;
    x1_norm=nan;
    x2_norm=nan;
    p1_norm=nan;
    p2_norm=nan;
    x_all=nan;
    p_all=nan;
    y2=nan;
    return
end
%raster sorted by dk
%sorted by touch order in every trial
[deltak_ver,idx]=sort(abs(deltak_hor));
raster_delta2=spikes(idx,:);
touch_idx2=touch_idx(idx);

if ~normalised
    x1_norm=0;
    x2_norm=0;
    p1_norm=0;
    p2_norm=0;
end

if ~do_plot
    
    [~,x1,p1]=plotTCurve(deltak_1,FR_1,'',N_first,percentile_1);
    [~,x2,p2,y2]=plotTCurve(deltak_2,FR_2,'',N_latter,percentile_2);
    [~,x_all,p_all]=plotTCurve(deltak_1,FR_1wong,'',N_first,percentile_1);
    if normalised
        [~,x1_norm,p1_norm]=plotTCurve(deltak_1_norm,FR_1,'',N_first,percentile_1);
        [~,x2_norm,p2_norm,y2_norm]=plotTCurve(deltak_2_norm,FR_2,'',N_latter,percentile_2);
    end
    
    


% subplot(3,3,6) 
%     hold on 
%     delta2=abs(deltak_hor(touch_idx==2));
%     delta3=abs(deltak_hor(touch_idx==3));
%     delta4=abs(deltak_hor(touch_idx>3));
%     errorbar([1 2 3 4],[mean(deltak_1) mean(delta2) mean(delta3) mean(delta4)],[std(deltak_1) std(delta2) std(delta3) std(delta4)])
%     
%     
 %    subplot(3,3,9) 
%     hold on 
%     FR_3=sum(spikes(touch_idx==3,100:130),2);
%     FR_4=sum(spikes(touch_idx>3,100:130),2);
%     errorbar([1 2 3 4],[mean(FR_1) mean(FR2) mean(FR_3) mean(FR_4)],[std(FR_1) std(FR2) std(FR_3) std(FR_4)])
%     
end
if do_plot
      
%     subplot(3,3,1)
%     plot(-100:100,mean(spikes(touch_idx==1,:)),'r')
%     hold on
%     plot(-100:100,mean(spikes(touch_idx>1,:)),'k')
%     xlabel('Time to touch onset [ms]')
%     legend('First touch','Latter touches')
%     ylabel('FR')
%     box off
    
%     subplot(3,3,2)
%     hold on
%    for i=1:size(deltak_ver,2)
%         s1=find(raster_delta2(i,100:end));
%         if ~isempty(s1)
%             if touch_idx2(i)==1
%                 plot(s1,i*ones(size(s1,2)),'.k'),hold on,
%             else
%                 plot(s1,i*ones(size(s1,2)),'.k'),hold on,
%             end
%         end
%         plot(sum(raster_delta2(i,100:130)),i,'o')
%     end
%         responses=sum(raster_delta2(:,100:130),2);
%      plot(responses,1:size(deltak_ver,2),'o')
%      hold on
%     xlabel('Time to touch onset [ms]')
%     h=histogram(1:size(deltak_ver,2),10,'Visible','off');
%     h.BinEdges(end)=size(deltak_ver,2);
%     h.BinEdges(1)=1;
%     for r=1:10
%         h.BinEdges(r)
%         r_mean(r)=mean(responses(h.BinEdges(r):h.BinEdges(r+1)));
%         x_rmean(r)=mean(h.BinEdges(r),h.BinEdges(r+1));
%     end
%     plot(r_mean,x_rmean)
    %xlim([0 100])
%     ax1 = gca; % current axes
%     ax1_pos = ax1.Position; % position of first axes
%     ax2 = axes('Position',ax1_pos,...
%         'XAxisLocation','top',...
%         'YAxisLocation','right',...
%         'Color','none','XColor','b','YColor','none');
%     hold on
%     zeropoint=0;
%     line(deltak_ver,1:size(deltak_ver,2),'Color','b','LineWidth',1.5)
%     line([zeropoint zeropoint],[0 size(deltak_ver,2)],'Color','b','LineWidth',1.5,'LineStyle','--')
%     line([0 size(deltak_ver,2)],[zeropoint zeropoint],'Color','b','LineWidth',1.5,'LineStyle','--')
%     xlim([min(deltak_ver) max(deltak_ver)])
%    title('Sorted by  \Delta \kappa');

%      line(1:size(deltak_ver,2),deltak_ver,'Color','b','LineWidth',1.5)
% %     line([zeropoint zeropoint],[0 size(deltak_ver,2)],'Color','b','LineWidth',1.5,'LineStyle','--')
%      line([0 size(deltak_ver,2)],[zeropoint zeropoint],'Color','b','LineWidth',1.5,'LineStyle','--')
%      ylim([min(deltak_ver) max(deltak_ver)])
%      title('Sorted by  \Delta \kappa');

    
    
    %subplot(3,3,3)
    %h1=histogram(deltak_1,10,'Normalization','probability','Facecolor','r');
    %histogram(deltak_2,'BinWidth',h1.BinWidth,'Normalization','probability','Facecolor','k')
%     delta2=abs(deltak_hor(touch_idx==2));
%     delta3=abs(deltak_hor(touch_idx==3));
%     delta4=abs(deltak_hor(touch_idx>3));
%     errorbar([1 2 3 4],[mean(deltak_1) mean(delta2) mean(delta3) mean(delta4)],[std(deltak_1) std(delta2) std(delta3) std(delta4)])
%     
%     box off
%     xlabel('\Delta \kappa')
%     ylabel('Normalised frequency')
%     plot(deltak_ver,1:size(deltak_ver,2),'k')
%     hold on 
    
    
    subplot(7,5,sub_plot_idx)
    hold on
    [nsamples2,x2,p2,y2]=plotTCurve(deltak_2,FR_2,'k',N_latter,percentile_2);    
    [~,x_all,p_all]=plotTCurve(deltak_1,FR_1wong,'',N_first,percentile_2);
    [nsamples1,x1,p1]=plotTCurve(deltak_1,FR_1,'r',N_first,percentile_1);
    [~,xwong,pwong]=plotTCurve(deltak_1,FR_1wong,'c',N_first,percentile_1);
    xlabel('\Delta \kappa')
    ylabel('FR')
    box off
    
%     subplot(3,3,5)
%     plot(nsamples1,'r')
%     hold on
%     plot(nsamples2,'k')
%     box off
%     xlabel('Segment')
%     ylabel('Nsamples')
%  
    
end


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
% means
P=polyfit(xtemp,y,1);
%all points
%P=polyfit(Amp,Nspikes,1);
f = polyval(P,xtemp);


if strcmp(plot_d,'')
else
    
    errorbar(xtemp,y,err,[ '.' plot_d])
    hold on
    plot(xtemp,f,['-o' plot_d])
    title(num2str(P(1)))
end
end