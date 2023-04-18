function sensitive_per_whisker
%% if neuron i is more sensitive to whisker 1, then create a tuning curve for first and later touches just for that whisker 
[~,~,~,~,~,~,deltak,touch_idx,FR,~,~,whisker]=Tcurve_touches_dk(Data.touch,Data.touch_per_whisker,total_psth,Data.deltak_w1,Data.deltak_w2,do_extra_plot,NpointsTcurve,0,percentile_1,percentile_2);
       
%% create tuning curve for first touch of the selected whisker 
pref_w=1; % or 2
selected_touch=touch_idx==1 & whisker==pref_w;
nsamples_1=sum(selected_touch);
[~,x1,p1]=plotTCurve(deltak(selected_touch),FR(selected_touch),'',N_first,percentile_1);

selected_touch=touch_idx>1 & whisker==pref_w;
nsamples_2=sum(selected_touch);
[~,x2,p2]=plotTCurve(deltak(selected_touch),FR(selected_touch),'',N_first,percentile_1);

end