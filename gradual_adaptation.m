function gradual_adaptation(fig4)
figure(fig4)
k3d=2;%1= k3d, 0=k_hor, 2=newk
counter=0;
av_unit=0;
folder=[{'Glu32_19092017H'};{'Glu32_21092017H'};{'Glu43_22122017H'};{'Glu35_10112017H'};{'Glu35_13112017H_1'};{'Glu35_13112017H_2'}];
%folder2=[{'Glu32_190917a_20170919_'},{'Glu32_210917a_20170921_'},{'Glu43_221217a_20171222_'},{'Glu35_101117a_20171110_'},{'Glu35_131117a_20171113_'},{'Glu35_131117a_20171113_'}];
[~,txt]=xlsread('Adaptation_units_list.xlsx',1);
%colours=[{'b'};{'g'};{'r'};{'c'};{'m'};{'k'}];
for f=1:size(folder,1)
    
    cd(folder{f})
    
    if k3d==1
        %k3D
        %         load('deltak_c2.mat')
        %         k_c2=touches_matrix/0.047;
        %         load('deltak_c1.mat')
        %         k_c1=touches_matrix/0.047;
        [~,~,~,matrix_norm]=normalised_by_trial(15,0);
        k_c1=matrix_norm(:,:,1);
        k_c2=matrix_norm(:,:,2);
    elseif k3d==0
        %k_h
        load('kinematicsgoc2.mat')
        k_c2=azimuth;
        load('kinematicsgoc1.mat')
        k_c1=azimuth;
        
    else

        load('newk1.mat')
        k_c1=newk1*1.7/0.047;
        load('newk2.mat')
        k_c2=newk2/0.047;
    end
    
    
    
    load('all_touches.mat')
    load('touches_whisker.mat')
    ff=dir(['*neural_data_S*']);
    clear total_psth
    total_psth=zeros(size(touches_matrix));
    for unit=1:size(ff,1)
        
        if onthelist(txt,[folder{f} ff(unit).name])
            load(ff(unit).name)
            total_psth=psth;
            
            
            delete=size(total_psth,2)-size(k_c1,2);
            k_c1=[zeros(size(k_c1,1),delete),k_c1];
            k_c2=[zeros(size(k_c1,1),delete),k_c2];
            
            clear deltak touch_idx FR
            
            %[x1,p1,x2,p2,x1_norm,p1_norm,x2_norm,p2_norm,x_all,p_all,deltak,touch_idx,FR,~,FR_prev,whisker]=Tcurve_touches_dk(touches_matrix,touches_whisker,total_psth,k_c1,k_c2,0,4);
            [x1,x2,x3,x4,p,deltak,touch_idx,FR,whisker]=Tcurve_touches_dk_per_touch(touches_matrix,touches_whisker,total_psth,k_c1,k_c2,4);
            [h_atte,pval_atte]=ttest2(FR(touch_idx==1),FR(touch_idx>1));
            if h_atte
                av_unit=av_unit+1;

                slope(av_unit,:)=p(:,1)';
                hold on
                ylabel('Slope')
                box off
                xlim([0.5 4.5])
                
                subplot(2,4,8)
                plot(p(:,2),'-','Color',[0.5 0.5 0.5])
                intercept(av_unit,:)=p(:,2)';
                hold on

                
                subplot(2,4,7)
                FR1_all=simspikes(deltak(touch_idx==1)',p(4,:),1); %predict first touch responses with later curves
                FR2_all=simspikes(deltak(touch_idx==2)',p(4,:),1); %predict second touch responses with later curves
                FR3_all=simspikes(deltak(touch_idx==3)',p(4,:),1); %predict second touch responses with later curves
                FR4_all=simspikes(deltak(touch_idx==4)',p(4,:),1); %predict second touch responses with later curves
                mean_pred=[mean(FR1_all(:)) mean(FR2_all(:)) mean(FR3_all(:)) mean(FR4_all(:))];
                mean_FR=[mean(FR(touch_idx==1)), mean(FR(touch_idx==2)),mean(FR(touch_idx==3)),mean(FR(touch_idx==4))];
                ratioFR(av_unit,:)=mean_pred./mean_FR;
                plot(mean_pred./mean_FR,'-','Color',[0.5 0.5 0.5])
                hold on
                clear FR p deltak touch_idx whisker
            end
        end
    end
    cd ..
    
end
subplot(2,4,7)
errorbar(1:4,mean(ratioFR),std(ratioFR),'.-k','LineWidth',2)
xlim([0.5 4.5])
xlabel('Touch number')
ylabel('Predicted FR/ Observed FR')
box off
xlim([0.5 4.5])

subplot(2,4,8)
errorbar(1:4,mean(intercept),std(intercept),'.-k','LineWidth',2)
xlabel('Touch number')
ylabel('Intercept')
xlim([0.5 4.5])
box off

% subplot(2,2,1)
% errorbar([1:4],mean(slope),std(slope),'.-k','LineWidth',2)

end

function FR=simspikes(x,p,Nsim)
lambda=p(1)*x+p(2);
lambda=repmat(lambda,1,Nsim);
FR = lambda;
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