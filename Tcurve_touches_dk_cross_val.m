function [lambda1,lambda2]=Tcurve_touches_dk_cross_val(touches,touches_whisker,psth,k_c1,k_c2,do_plot,N_first)
%this code assumes that the input is dk1 and dk2
counter=1;
k_hor(:,:,1)=k_c1;
k_hor(:,:,2)=k_c2;
window=100;
window_ad=600;
normalised=1;
%N_first=5;
N_latter=N_first;
percentile_1=1;
percentile_2=0.95;
if do_plot
    figure
end
    
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
                deltak_hor(counter)=abs(mean(k_hor(trial,idx(t):idx(t)+5,whisker(counter))));%%curvature hor   
            else
                deltak_hor(counter)=abs(mean(k_hor(trial,idx(t):end,whisker(counter))));%%curvature hor
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
FR=sum(spikes(:,100:100+upto),2);
deltak=abs(deltak_hor);
deltak_norm=abs(deltak_hor)./std_dk;

%first touch
FR_1=sum(spikes(touch_idx==1,100:100+upto),2);
deltak_1=abs(deltak_hor(touch_idx==1));
%for the first touch, std of delta kappa should be zero by definition

%latter touches
FR_2=sum(spikes(touch_idx>1,100:100+upto),2);
deltak_2=abs(deltak_hor(touch_idx>1));

FR2=sum(spikes(touch_idx==2,100:100+upto),2);
total_sp=sum(FR_1)+sum(FR2);
if total_sp<20
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
    lambda1=nan;
    lambda2=nan;
    return
end




if ~do_plot
    %% cross-validation for first touch
    for i=1:numel(deltak_1) 
        cross_idx1=1:numel(deltak_1);
        cross_idx1(i)=[];
    [~,~,p1(i,:)]=plotTCurve(deltak_1(cross_idx1),FR_1(cross_idx1),'',N_first,percentile_1);
    %%prediction for the excluded sample
    lambda1(i)=p1(i,1)*deltak_1(i)+p1(i,2);
    end
    
    %%% cross_validaiton for first touch
    for i=1:numel(deltak_2) 
    cross_idx2=1:numel(deltak_2);
        cross_idx2(i)=[];
    [~,~,p2(i,:)]=plotTCurve(deltak_2(cross_idx2),FR_2(cross_idx2),'',N_latter,percentile_2);
    lambda2(i)=p2(i,1)*deltak_2(i)+p2(i,2);
    end
end

if do_plot
 
    
    
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