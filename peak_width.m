function [peakW,peak_latency]=peak_width(spikes,touch_idx)

peakW=zeros(1,4);
peak_latency=zeros(1,4);

boxcar=[zeros(1,5) ones(1,5)];

for itouch=1:4

    if itouch==4
        filteredFR=conv(mean(spikes(touch_idx>=itouch,100:end)),boxcar,'same');
    else
        filteredFR=conv(mean(spikes(touch_idx==itouch,100:end)),boxcar,'same');
    end

    [MaxF,idx_max]=max(filteredFR(1:end));
    end_peak=find(filteredFR(idx_max:end)<=MaxF/2.5,1,'first');
    if isempty(end_peak)
        end_peak=nan;
    end
    start_peak=idx_max-find(filteredFR(1:idx_max)<=MaxF/2.5,1,'last');

    if isempty(start_peak)
        start_peak=nan;
    end
    
    subplot(2,3,6)
    hold on
    plot(filteredFR)

    peakW(itouch)=end_peak+start_peak;
    peak_latency(itouch)=idx_max;

end

end