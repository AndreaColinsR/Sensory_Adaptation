function first_FR=cross_whisker_adapt_v2
fig1=figure;
fig2=figure;
av_unit=1;
folder=[{'Glu43_22122017H'};{'Glu32_19092017H'};{'Glu32_21092017H'};{'Glu35_10112017H'};{'Glu35_13112017H/first_segment'};{'Glu35_13112017H/second_segment'}];
[~,txt]=xlsread('./Adaptation/Adaptation_units_list.xlsx',1);
for f=1:size(folder,1)-1
    cd([folder{f} '/modelling/'])
    
    
    window=100;
    ff=dir('*neural_data_S*');
    load(ff(1).name,'psth')
    total=zeros(size(psth));
    load('newk1.mat','newk1')
    load('newk2.mat','newk2')
    newk=newk1*1.7/0.047;
    newk(:,:,2)=newk2/0.047;
    for unit=1:size(ff,1)
        if onthelist(txt,[folder{f} ff(unit).name]) 
            load(ff(unit).name,'psth')
            %total=total+psth;
            %end
            %psth=total;
            load('all_touches.mat','touches_matrix')
            load('touches_whisker.mat','touches_whisker')
           
            counter=1;
            counter2=1;
            final=size(psth,1);
            if contains(pwd,'Glu43_22122017H')
                %pwd
                final=83;
            end
            
            if contains(pwd,'Glu32_21092017H')
                %pwd
                final=97;
            end
            
            for trial=1:final
                idx1=find(diff(touches_whisker(trial,:,1))>0);
                idx2=find(diff(touches_whisker(trial,:,2))>0);
                [idx,ind]=sort([idx1 idx2]);
                idx=idx+1;
                whisker_trial=[ones(size(idx1)) ones(size(idx2))*2];
                whisker_trial=whisker_trial(ind);
                %          plot(touches_whisker(trial,:,1),'b')
                %          hold on
                %          plot(touches_whisker(trial,:,2),'g')
                %          plot(total(trial,:),'k')
                %          hold off
                %select for appropiate trials type
                
                if size(idx,2)>1 && (whisker_trial(1)==whisker_trial(2)) && idx(2)-idx(1)>20
                    for t=1:2
                        [end2,a]=min([idx(t)+window size(psth,2)]);
                        
                        dk(counter)=abs(mean(newk(trial,idx(t):idx(t)+5,whisker_trial(t))));
                        
                        
                        if a==1
                            spikes(counter,:)=psth(trial,idx(t)-100:idx(t)+window);
                            
                        else
                            aux=zeros(1,size(spikes(1,:),2)-size(psth(trial,idx(t)-100:end2),2));
                            spikes(counter,:)=[psth(trial,idx(t)-100:end2) aux];
                            
                        end
                        
                        touch_idx(counter)=t;
                        counter=counter+1;
                    end
                    
                    t=find(whisker_trial(1)-whisker_trial,1,'First');
                    %exclude two whiskers touching at the same time
                    
                    if ~isempty(t) && (sum(touches_whisker(trial,idx(t),:))<2)
                        [end2,a]=min([idx(t)+window size(psth,2)]);
                        dk2(counter2)=abs(mean(newk(trial,idx(t):idx(t)+5,whisker_trial(t))));
                        whisker_test(counter2)=whisker_trial(t);
                        %             subplot(2,1,1)
                        %             plot([dk(counter-2) dk(counter-1) dk2(counter2)])
                        
                        if a==1
                            spikes2(counter2,:)=psth(trial,idx(t)-100:idx(t)+window);
                            
                        else
                            aux=zeros(1,size(spikes2(1,:),2)-size(psth(trial,idx(t)-100:end2),2));
                            spikes2(counter2,:)=[psth(trial,idx(t)-100:end2) aux];
                            
                        end
                        
                        %             subplot(2,1,2)
                        %             plot([ sum(spikes(counter-2,100:130)) sum(spikes(counter-1,100:130)) sum(spikes2(counter2,100:130))])
                        %             trial
                        %             pause
                        counter2=counter2+1;
                    end
                end
                clear idx whisker_trial
            end
            
            FR1=sum(spikes(touch_idx==1,100:130)');
            FR2=sum(spikes(touch_idx==2,100:130)');
            FR3=sum(spikes(touch_idx==3,100:130)');
            FR1_1=sum(spikes2(:,100:130)');
           
 
            deltak1=dk(touch_idx==1)';
            deltak2=dk(touch_idx==2)';
            %deltak3=dk(touch_idx==3)';
            
            
            
            
            [~,p1,~,p2,~,~,~,~,~,~,deltaktotal,touch_idxtotal,FRtotal,~,~,whisker]=Tcurve_touches_dk(touches_matrix,touches_whisker,psth,newk(:,:,1),newk(:,:,2),0,4);
            f1_first_touch = polyval(p1,dk2);

            deltaktotal(touch_idxtotal<3)=[];
            FRtotal(touch_idxtotal<3)=[];
            whisker(touch_idxtotal<3)=[];
            touch_idxtotal(touch_idxtotal<2)=[];

            for n_ex=1:numel(dk2)
                idx_tmp(n_ex)=find(deltaktotal==dk2(n_ex));
            end
            
            deltaktotal(idx_tmp)=[];
            FRtotal(idx_tmp)=[];
            touch_idxtotal(idx_tmp)=[];
            whisker(idx_tmp)=[];
            
            
            %produce one tuning curve per whisker
            pw1=0;
            pw2=0;
            fw1=[];
            fw2=[];

            if any(whisker_test==1)

            [~,~,pw1]=plotTCurve(deltaktotal(whisker==1),FRtotal(whisker==1),'',4,0.95);
            fw1 = polyval(pw1,dk2(whisker_test==1));
            end

            if any(whisker_test==2)
            [~,~,pw2]=plotTCurve(deltaktotal(whisker==2),FRtotal(whisker==2),'',4,0.95);
            fw2 = polyval(pw2,dk2(whisker_test==2));
            end
            %the prediction from the tuning curves for N>=3
            f1=[fw1,fw2]';


            figure(fig1)
            subplot(3,2,[1 2])
            %errorbar(1:4,[mean([FR1',FR2',FR3']) mean(FR1_1)] ,[std([FR1',FR2',FR3']) std(FR1_1)],'.-')
            FR_touch(av_unit,:)=[mean([FR1',FR2']) mean(FR1_1)];
            
            [h(av_unit),p(av_unit)] = ttest2(FR2',FR1_1,'Tail','left');
            if h(av_unit)>0
              errorbar(1:3,[mean([FR1',FR2']) mean(FR1_1)],[std([FR1',FR2']) std(FR1_1)],'b')
              
            else
              errorbar(1:3,[mean([FR1',FR2']) mean(FR1_1)],[std([FR1',FR2']) std(FR1_1)],'k')
           
            end
            hold on
            xticks([1 2 3 4])
            xticklabels({'1st touch','2nd touch','3rd touch','1st touch'})
            box off
            xlim([ 0.8 4.2])
            ylabel('Firing rate')
            
            
           
            
            first_FR(av_unit,:)=[mean(f1), mean(FR1_1)];
            
            [h_pred1(av_unit),p_pred1(av_unit)] = ttest2(FR1_1,f1,'Tail','right');

            musall_style(av_unit)=(mean(FR1_1)-mean(f1))/(mean(FR1_1)+mean(f1))
            if h_pred1(av_unit)>0
                 subplot(3,2,5)
                 hold on 
                plot([1 2],first_FR(av_unit,:),'-og')
                
                subplot(3,2,6)
                 hold on 
                 h_g=ttest2(f1_first_touch,FR1_1);
                 if h_g>0
                  plot([1 2],[mean(f1_first_touch) first_FR(av_unit,2)],'-og')
                 else
                  plot([1 2],[mean(f1_first_touch) first_FR(av_unit,2)],'-oc')   
                 end
                  
            else
                subplot(3,2,5)
                hold on
                 plot([1 2],first_FR(av_unit,:),'-oy')
            end

            %plot(mean(f1), mean(FR1_1),'-ob')
            title('Prediction first touch')
            hold on
            %plot([1 2],[mean(f1_control),mean(FR1)],'-ok')
            %plot(mean(f1_control),mean(FR1),'-ok')
            box off
            xlabel('Prediction tuning curve')
            ylabel('real FR')
            
            subplot(3,2,6)
%             second_FR(av_unit,:)=[mean(f2) mean(FR1_1)];
% 
%             [h_pred(av_unit),p_pred(av_unit)] = ttest2(FR1_1,f2,'Tail','right');
%             if h_pred(av_unit)>0
%             plot([1 2],second_FR(av_unit,:),'-og')
%             else
%             plot([1 2],second_FR(av_unit,:),'-or')    
%             end
%             
            hold on
            %tmp(av_unit,:)=[mean(f2_control), mean(FRtotal(touch_idxtotal>=3))];
            %plot([1 2],[mean(f2_control), mean(FRtotal(touch_idxtotal>=3))],'-ok')

            box off
            xlabel('Prediction tuning curve')
            ylabel('real FR')
            title('Prediction later touches')
            av_unit=av_unit+1;
        end
        
        
        clear f1_first_touch idx_tmp FR_1 FR1 FR2 FR3  spikes spikes2 psth total touch_idx separated_FR together_FR separated_dk together_dk
        
    end
    
    deltak_touch(f,:)=[mean([deltak1,deltak2]) mean(dk2)];
    
    samples(f)=numel(dk2);
    %[hd1,pd1] = ttest2(deltak1,dk2,'Tail','right');
    %[hd2,pd2] = ttest2(deltak2,dk2);
    
    subplot(3,2,[3 4])
    errorbar(1:3,deltak_touch(f,:),[std([deltak1,deltak2]) std(dk2)],'b')
    hold on
    xticks([1 2 3 4])
    xticklabels({'1st touch','2nd touch','3rd touch','1st touch'})
    box off
    xlim([ 0.8 4.2])
    ylabel('\kappa [1/mm]')
    clear deltak1 deltak2  deltak3 dk2 whisker_test
    
%     subplot(3,2,[1 2])
%     hold off
%     pause
    
%     figure(fig2)
%     if f==2
%         [separated_FR,together_FR,separated_dk,together_dk]=integration_whisker_t(1);
%         
%     end
%     subplot(3,3,7)
%     %plot([1 2],[mean(separated_FR)/mean(separated_FR) mean(together_FR)/mean(separated_FR)],'-o')
%     hold on
%     box off
%     xlim([0.8 2.2])
%     xticks([1 2])
%     %ylim([0 2])
%     xticklabels({'\Delta t> 30 ms','\Delta t< 20 ms'})
%     ylabel('FR/FR separated')
%     title('Time difference between touches')
    
%     subplot(3,3,8)
%     plot([1 2],[mean(sum(separated_dk,2))/mean(sum(separated_dk,2)) mean(sum(together_dk,2))/mean(sum(separated_dk,2))],'-o')
%     hold on
%     box off
%     xlim([0.8 2.2])
%     xticks([1 2])
%     ylim([0 2])
%     xticklabels({'\Delta t> 30 ms','\Delta t< 20 ms'})
%     ylabel('FR/FR prediction')
    
    
    %predict fr with dk
%     f1_1 = polyval(p2,separated_dk(:,1));
%     f1_2 = polyval(p2,separated_dk(:,2));
%     pred_separated=f1_1+f1_2;
%     
%     f2_1 = polyval(p2,together_dk(:,1));
%     f2_2 = polyval(p2,together_dk(:,2));
%     pred_together=f2_1+f2_2;
%     
%     subplot(3,3,9)
%     plot([1 2],[mean(separated_FR)/mean(pred_separated) mean(together_FR)/mean(pred_together)],'-o')
%     hold on
%     box off
%     xlim([0.8 2.2])
%     xticks([1 2])
%     ylim([0 2])
%     xticklabels({'\Delta t> 30 ms','\Delta t< 20 ms'})
%     ylabel('FR/FR prediction')
    
    if f==size(folder,1)-1
        cd ..
        cd ..
        cd ..
    else
        cd ..
        cd ..
        
    end
   
end
[h_pred1 (av_unit-1)]
positive_fraction=sum(h_pred1)/(av_unit-1)
musall_style
disp(['samples=' num2str(mean(samples)) '+-' num2str(std(samples))])
mean_p_all_units=mean(p_pred1)
mean_p_positive=mean(p_pred1(h_pred1>0))

figure(fig1)

%subplot(3,2,[1 2])
%errorbar(1:3,mean(FR_touch),std(FR_touch),'.-k')



subplot(3,2,[3 4])
%errorbar(1:3,mean(deltak_touch),std(deltak_touch),'.-k')

subplot(3,2,6)
%plot([0 60],[0 60])
%errorbar([1 2],mean(second_FR),std(second_FR),'.-c')
%xlim([0.8 2.2])
%errorbar(mean(second_FR(:,1)),mean(second_FR(:,2)),std(second_FR(:,1)),'horizontal','.b')
%errorbar(mean(second_FR(:,1)),mean(second_FR(:,2)),std(second_FR(:,2)),'vertical','.b')

% scatterDiagHist(second_FR(:,1),second_FR(:,2))
% hold on
% scatterDiagHist(tmp(:,1),tmp(:,2))



subplot(3,2,5)
%plot([0 60],[0 60])
errorbar([1 2],mean(first_FR),std(first_FR),'.-m')
xlim([0.8 2.2])
%errorbar(mean(first_FR(:,1)),mean(first_FR(:,2)),std(first_FR(:,1)),'horizontal','.b')
%errorbar(mean(first_FR(:,1)),mean(first_FR(:,2)),std(first_FR(:,2)),'vertical','.b')
%dim = [.7 0.1 .1 .1];
subplot(3,2,6)
xlim([0.8 2.2])
%figure
%%real FR

% subplot(3,1,1)
% [~,idx_FR]=sort(FR_touch(:,3)./FR_touch(:,2));
% 
% plot(FR_touch(idx_FR,3)./FR_touch(idx_FR,2),'o')
% hold on
% plot([0 size(FR_touch,1)],[1 1],'k')
% xlabel('Units')
% ylabel('FR First touch new w/FR Second touch')

% subplot(3,1,2)
% sum(h)
% %sum(h_pred & h)
% plot(p(idx_FR),'o')
% hold on
% plot([0 size(FR_touch,1)],[0.05 0.05],'k')
% xlabel('Units')
% ylabel('P-value equal mean')

% 
% %%population level
% subplot(3,1,3)
% histogram(FR_touch(:,2),1:10:80)
% hold on
% histogram(FR_touch(:,3),1:10:80)
% [h,p] = ttest2(FR_touch(:,2),FR_touch(:,3),'Tail','left');
% str=['p-value= ' num2str(p,3)];
% annotation('textbox',dim,'String',str,'FitBoxToText','on');
% xlabel('Mean firing rate [Hz] per unit')
% ylabel('Number of units')
% legend('Second Touch','First touch new w')
% 
% 
% [hd,pd] = ttest2(deltak_touch(:,2),deltak_touch(:,3),'Tail','right');
% 
% figure
% %%real FR and predicted FR
% %%population level
% subplot(3,1,3)
% %histogram(second_FR(:,1),1:10:80)
% hold on
% histogram(FR_touch(:,3),1:10:80)
% 
figure
musall_style(musall_style==-1)=[];
histogram(musall_style,-1:0.1:1)
sum(musall_style>0)/numel(musall_style)
[h_musall,p_musall]=ttest(musall_style)
p = signrank(musall_style)
end

function [separated_FR,together_FR,separated_dk,together_dk]=integration_whisker_t(do_plot)
window=530;
ff=dir('*neural_data_S*');
load(ff(1).name,'psth')
total=zeros(size(psth));
load('newk1.mat','newk1')
load('newk2.mat','newk2')
newk=newk1*1.7/0.047;
newk(:,:,2)=newk2/0.047;

for unit=1:size(ff,1)
    load(ff(unit).name,'psth')
    total=total+psth;
end
psth=total;
load('all_touches.mat')
load('touches_whisker.mat')
touches=touches_matrix;
counter=1;
for trial=1:size(total,1)
    idx1=find(diff(touches_whisker(trial,:,1))>0);
    idx2=find(diff(touches_whisker(trial,:,2))>0);
    whisker_t=ones(size(idx1));
    whisker_t2=ones(size(idx2))*2;
    [idx,ind]=sort([idx1 idx2]);
    whisker_trial=[ones(size(idx1)) ones(size(idx2))*2];
    whisker_trial=whisker_trial(ind);
    
    %select for appropiate trials type
    for t=3:size(idx,2)

        if (whisker_trial(t-1)~=whisker_trial(t)) && idx(t)-idx(t-1)<=500 && idx(t)<3483 && (sum(touches_whisker(trial,idx(t)+1,:))<2)
            %detach(counter)=find(diff(touches_whisker(trial,idx(1):end,whisker_trial))<0,1,'First');
            [end2,a]=min([idx(t-1)+window size(psth,2)]);
            dk1(counter)=mean(abs(newk(trial,idx(t-1):idx(t-1)+5,whisker_trial(t-1))));
            
            dk2(counter)=mean(abs(newk(trial,idx(t):idx(t)+5,whisker_trial(t))));
            if a==1
                spikes(counter,:)=psth(trial,idx(t-1)-100:idx(t-1)+window);
                
            else
                aux=zeros(1,size(spikes,2)-size(psth(trial,idx(t-1)-100:end2),2));
                spikes(counter,:)=[psth(trial,idx(t-1)-100:end2) aux];
                
            end
            second_t(counter)=idx(t)-idx(t-1);
            counter=counter+1;
            
        end
    end
    
end

%sort spikes according to intertouch interval
[second_t,idx_s]=sort(second_t);
spikes=spikes(idx_s,:);
dk1=dk1(idx_s);
dk2=dk2(idx_s);
%detach=detach(idx_s);
%make rasterplot
if do_plot
    subplot(3,3,[1:6])
    tmp=find(second_t>30,1,'First');
    for i=tmp:size(spikes,1)
        t=find(spikes(i,:));
        for j=1:numel(t)
            plot(t(j)-100,i,'.k')
            hold on
        end
        
        plot([second_t(i) second_t(i) ],[i-0.5 i+0.5],'r')
        plot([second_t(i)+30 second_t(i)+30],[i-0.5 i+0.5],'r')
        %plot([detach(i) detach(i) ],[i-0.5 i+0.5],'g')
        %plot([30 30],[i-0.5 i+0.5],'g')
    end
    title('Same whisker first and second touch')
    plot([0 0],[tmp size(spikes,1)],'g')
    plot([30 30],[tmp size(spikes,1)],'g')
    xlabel('Time from first touch')
    box off
    ylabel('trial number')
    ylim([0 size(spikes,1)+0.1])
end

%compare FRs
subplot(3,3,[7])
%separated first touch
from=find(second_t>30,1,'First');
FR1=sum(spikes(from:end,100:130),2);
FR3=[];
for i=from:numel(second_t)
    FR3(i-from+1,1)=sum(spikes(i,second_t(i)+100:second_t(i)+130),2);
end

errorbar([1 2],[mean(FR1) mean(FR3)],[std(FR1) std(FR3)])
%mixed responses 1+2 touch deltat<15
to=find(second_t<=20,1,'Last');
FR2=sum(spikes(1:to,100:150),2);
%errorbar([1 2],[mean(FR1+FR3) mean(FR2)],[std(FR1+FR3) std(FR2)],'-.k')
separated_FR=FR1+FR3;
separated_dk=[dk1(from:end)' dk2(from:end)'];
together_FR=FR2;
together_dk=[dk1(1:to)' dk2(1:to)'];

% box off
% hold on
%
%
% xticks([ 0 1 2 3 4 5])
% xticklabels({'','First touch','\Delta t<15','30<\Delta t <70','70<\Delta t <150','150<\Delta t'})
% ylabel('Mean Firing rate [sp/stimulus]')

%%dk
% subplot(3,3,[3 6])
% plot(dk1,1:size(spikes,1),'g')
% hold on
% plot(dk2,1:size(spikes,1),'r')
% box off
% xlabel('dk')

end

function present=onthelist(txt,unit)

for i=1:size(txt,1)
    Index=strcmp([txt{i,:}],unit);
    if Index==0
        present=0;
    else
        present=1;
        return
    end
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