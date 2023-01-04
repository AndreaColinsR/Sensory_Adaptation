function cross_whisker_adaptation(fig5)

av_unit=1;
folder=[{'Glu43_22122017H'};{'Glu32_19092017H'};{'Glu32_21092017H'};{'Glu35_10112017H'};{'Glu35_13112017H_1'};{'Glu35_13112017H_2'}];
[~,txt]=xlsread('Adaptation_units_list.xlsx',1);
percentile_1=0.95;
percentile_2=0.95;
window=100;
nunits=size(txt,1);
first_FR=nan(nunits,2);
h_pred1=nan(nunits,1);
p_pred1=nan(nunits,1);
WSAI=nan(nunits,1);

%% Andrea
% check why we are ignoring the last session
subplot(1,2,1)
hold on
for f=1:size(folder,1)
    cd([folder{f}])
    
    

    ff=dir('*neural_data_S*');
    load(ff(1).name,'psth')
    load('newk1.mat','newk1')
    load('newk2.mat','newk2')
    newk=newk1*1.7/0.047;
    newk(:,:,2)=newk2/0.047;
    
    for unit=1:size(ff,1)
        if onthelist(txt,[folder{f} ff(unit).name])
            load(ff(unit).name,'psth')
            load('all_touches.mat','touches_matrix')
            load('touches_whisker.mat','touches_whisker')
            
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
            
            %% Calculate test response
            for trial=1:final
                idx1=find(diff(touches_whisker(trial,:,1))>0);
                idx2=find(diff(touches_whisker(trial,:,2))>0);
                [idx,ind]=sort([idx1 idx2]);
                idx=idx+1;
                whisker_trial=[ones(size(idx1)) ones(size(idx2))*2];
                whisker_trial=whisker_trial(ind);
                %select for appropiate trials type
                
                if size(idx,2)>1 && (whisker_trial(1)==whisker_trial(2)) && idx(2)-idx(1)>20
                    
                    t=find(whisker_trial(1)-whisker_trial,1,'First');
                    %exclude two whiskers touching at the same time
                    
                    if ~isempty(t) && (sum(touches_whisker(trial,idx(t),:))<2)
                        [end2,a]=min([idx(t)+window size(psth,2)]);
                        dk2(counter2)=abs(mean(newk(trial,idx(t):idx(t)+5,whisker_trial(t))));
                        whisker_test(counter2)=whisker_trial(t);
                        
                        
                        if a==1
                            spikes2(counter2,:)=psth(trial,idx(t)-100:idx(t)+window);
                            
                        else
                            aux=zeros(1,size(spikes2(1,:),2)-size(psth(trial,idx(t)-100:end2),2));
                            spikes2(counter2,:)=[psth(trial,idx(t)-100:end2) aux];
                            
                        end
                        
                        counter2=counter2+1;
                    end
                end
                clear idx whisker_trial
            end
            
            % Test response
            FR1_1=sum(spikes2(:,100:130),2);
            
            %% Calculate prediction 
            
            [~,~,~,~,~,~,~,~,~,~,deltaktotal,touch_idxtotal,FRtotal,~,~,whisker]=Tcurve_touches_dk(touches_matrix,touches_whisker,psth,newk(:,:,1),newk(:,:,2),0,4,2,percentile_1,percentile_2);
            
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
            %% the prediction from the tuning curves for N>=3
            f1=[fw1,fw2]';

            
            %% Compare prediction to test response
            first_FR(av_unit,:)=[mean(f1), mean(FR1_1)];
            
            [h_pred1(av_unit),p_pred1(av_unit)] = ttest2(FR1_1,f1,'Tail','right');
            
            %% Compute Whisker Specific Adaptation Index  as in Musall et al. 
            
            WSAI(av_unit)=(mean(FR1_1)-mean(f1))/(mean(FR1_1)+mean(f1));
            
            
            if h_pred1(av_unit)>0
                
                plot([1 2],first_FR(av_unit,:),'-og')
                
            else
                plot([1 2],first_FR(av_unit,:),'-ok')
            end
            
            
            av_unit=av_unit+1;
        end
        
        
        clear idx_tmp FR_1  spikes spikes2 psth touch_idx touches_matrix touches_whisker
    end
    
    
    samples(f)=numel(dk2);
    
    
    clear dk2 whisker_test
    
    
    cd ..
    
    
end
positive_fraction=sum(h_pred1)/(av_unit-1);
disp(['samples=' num2str(mean(samples)) '+-' num2str(std(samples))])
mean_p_all_units=mean(p_pred1)
mean_p_positive=mean(p_pred1(h_pred1>0))

figure(fig5)
subplot(1,2,1)
errorbar([1 2],mean(first_FR,'omitnan'),std(first_FR,'omitnan'),'.-m')
xlim([0.8 2.2])
box off
ylabel('Firing Rate [spikes/touch]')
xticks([1 2])
xticklabels({'Predicted unadapted state','Test response'})


subplot(1,2,2)
WSAI(WSAI==-1)=[];
histogram(WSAI,-1:0.1:1)
hold on
plot([0 0],[0 4],'Color',[0.5 0.5 0.5])
box off
xlim([-1.1 1.1])
ylim([0 4])
xlabel('Whisker-Specific Adaptation Index')
ylabel('Number of units')

%% stats
%sum(musall_style>0)/numel(musall_style)
[h_musall,p_musall]=ttest(WSAI)
p = signrank(WSAI)

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