function cross_whisker_adaptation(fig5)
av_unit=1;
ff=dir('*.mat*');
Nsessions=size(ff,1);


percentile_FR=0.95;
window=100;
nunits=33;

first_FR=nan(nunits,2);
h_pred1=nan(nunits,1);
p_pred1=nan(nunits,1);
WSAI=nan(nunits,1);
samples=zeros(Nsessions,1);


%% Define population FR per mouse
mouse_id=zeros(Nsessions,1);

for f=1:Nsessions
    
mouse_id(f)=str2double(ff(f).name(4:5));
end

mouse_number=unique(mouse_id);
WSAI_mouse=cell(numel(mouse_number),1);



subplot(1,2,1)
hold on
for f=1:Nsessions
    
    
    load(ff(f).name,'Data')
    newk=Data.deltak_w1;
    newk(:,:,2)=Data.deltak_w2;
    
    this_mouse=str2double(ff(f).name(4:5))==mouse_number;
    
    for i_unit=1:size(Data.unit,2)
        
        counter2=1;
        final=size(Data.unit(i_unit).spikes,1);
        if strcmp(ff(f).name,'Glu43_22122017H.mat')
            final=83;
        end
        
        if strcmp(ff(f).name,'Glu32_21092017H.mat')
            
            final=97;
        end
        
        %% Calculate test response
        for trial=1:final
            idx1=find(diff(Data.touch_per_whisker(trial,:,1))>0);
            idx2=find(diff(Data.touch_per_whisker(trial,:,2))>0);
            [idx,ind]=sort([idx1 idx2]);
            idx=idx+1;
            whisker_trial=[ones(size(idx1)) ones(size(idx2))*2];
            whisker_trial=whisker_trial(ind);
            %select for appropiate trials type
            
            if size(idx,2)>1 && (whisker_trial(1)==whisker_trial(2)) && idx(2)-idx(1)>20
                
                t=find(whisker_trial(1)-whisker_trial,1,'First');
                %exclude two whiskers touching at the same time
                
                if ~isempty(t) && (sum(Data.touch_per_whisker(trial,idx(t),:))<2)
                    [end2,a]=min([idx(t)+window size(Data.unit(i_unit).spikes,2)]);
                    dk2(counter2)=abs(mean(newk(trial,idx(t):idx(t)+5,whisker_trial(t))));
                    whisker_test(counter2)=whisker_trial(t);
                    
                    
                    if a==1
                        spikes2(counter2,:)=Data.unit(i_unit).spikes(trial,idx(t)-100:idx(t)+window);
                        
                    else
                        aux=zeros(1,size(spikes2(1,:),2)-size(Data.unit(i_unit).spikes(trial,idx(t)-100:end2),2));
                        spikes2(counter2,:)=[Data.unit(i_unit).spikes(trial,idx(t)-100:end2) aux];
                        
                    end
                    
                    counter2=counter2+1;
                end
            end
            clear idx whisker_trial
        end
        
        % Test response
        FR1_1=sum(spikes2(:,100:130),2);
        
        %% Calculate prediction
        
        [~,~,~,~,~,~,deltaktotal,touch_idxtotal,FRtotal,~,~,whisker]=Tcurve_touches_dk(Data.touch,Data.touch_per_whisker,Data.unit(i_unit).spikes,newk(:,:,1),newk(:,:,2),0,4,2,percentile_FR,percentile_FR);
        deltaktotal=abs(deltaktotal);
        
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
            
            [~,~,pw1]=plotTCurve(deltaktotal(whisker==1),FRtotal(whisker==1),'',4,percentile_FR);
            fw1 = polyval(pw1,dk2(whisker_test==1));
        end
        
        if any(whisker_test==2)
            [~,~,pw2]=plotTCurve(deltaktotal(whisker==2),FRtotal(whisker==2),'',4,percentile_FR);
            fw2 = polyval(pw2,dk2(whisker_test==2));
        end
        %% the prediction from the tuning curves for N>=3
        f1=[fw1,fw2]';
        
        
        %% Compare prediction to test response
        first_FR(av_unit,:)=[mean(f1), mean(FR1_1)];
        
        [h_pred1(av_unit),p_pred1(av_unit)] = ttest2(FR1_1,f1,'Tail','right');
        
        %% Compute Whisker Specific Adaptation Index  as in Musall et al.
        
        WSAI(av_unit)=(mean(FR1_1)-mean(f1))/(mean(FR1_1)+mean(f1));
        WSAI_mouse{this_mouse}= [WSAI_mouse{this_mouse};WSAI(av_unit)];
        
        if h_pred1(av_unit)>0
            
            plot([1 2],first_FR(av_unit,:),'-og')
            
        else
            plot([1 2],first_FR(av_unit,:),'-ok')
        end
        
        
        av_unit=av_unit+1;
        
        
        clear idx_tmp FR_1  spikes spikes2 touch_idx
    end
    
    
    samples(f)=numel(dk2);
    
    
    clear dk2 whisker_test
    
    
end


disp(['samples=' num2str(mean(samples)) '+-' num2str(std(samples))])


%% test that FR is higher than predicted across neurons
[~,p_total]=ttest2(first_FR(:,1),first_FR(:,2));
% remove outlier and test again
[~,idx_max]=max(first_FR(:,2));
first_FR_outlier=first_FR;
first_FR_outlier(idx_max,:)=[];
[~,p_total_outlier]=ttest2(first_FR_outlier(:,1),first_FR_outlier(:,2));

figure(fig5)

subplot(1,2,1)
errorbar([1 2],mean(first_FR,'omitnan'),std(first_FR,'omitnan'),'.-m')
xlim([0.8 2.2])
box off
ylabel('Firing Rate [spikes/touch]')
xticks([1 2])
xticklabels({'Predicted unadapted state','Test response'})
text(1, 3, ['p-value (all) =' num2str(p_total,'%.2e')])
text(1, 2.5, ['p-value (w/o outlier) =' num2str(p_total_outlier,'%.2e')])

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
[~,p_musall]=ttest(WSAI);
p = signrank(WSAI);
text(-1, 3.5, ['p-value (signrank)= ' num2str(p)])
text(-1, 3, ['p-value (ttest)= ' num2str(p_musall)])
positive_I_ratio=sum(WSAI>0)/numel(WSAI)

disp('')
disp('')
disp('-----------------------------------------')
disp('Median I  per mouse')

for im=1:numel(mouse_number)
    
disp(['Mouse ' num2str(im) ' =' num2str(median(WSAI_mouse{im}))])

end
end