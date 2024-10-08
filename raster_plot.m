function raster_plot(spikes,touch_idx,nplot)
%check that both variables have the same number of rows
if size(spikes,1)~=size(touch_idx,2)
    display(['The variables have a different number of rows Spikes = ' num2str(size(spikes,1)) ' touch = ' num2str(size(touch_idx,2))])
    keyboard
end
%% plot rasterplot chronological order

subplot(4,3,nplot(1))
hold on
for i=1:size(spikes,1)
    s1=find(spikes(i,:))-100;
    if ~isempty(s1)
        plot(s1,i*ones(size(s1,2)),'.k'),
    end
end

plot([0 0],[0 size(spikes,1)],'Color',[0.5 0.5 0.5])
box off
xlabel('Time from touch onset [ms]')
ylabel('Touch number')
xlim([-100 100])
ylim([0 size(spikes,1)])
title('Chronological order')
hold off

% Filter the mean of the spikes
xtime=(1:size(spikes,2))-100;
boxcar=[zeros(1,5) ones(1,5)/5];
filteredFR=conv(mean(spikes,1),boxcar,'same');

subplot(4,3,nplot(3))
hold on
plot(xtime,filteredFR,'k')
    
plot([0 0],[0 2],'Color',[0.5 0.5 0.5])
box off
xlabel('Time from touch onset [ms]')
xlim([-100 100])
ylim([0 2])
hold off

%% Order within trials
[touch_idx,idx]=sort(touch_idx);
spikes=spikes(idx,:);

%% plot rasterplot by touch within trial
subplot(4,3,nplot(2))
hold on
for i=1:size(spikes,1)
    s1=find(spikes(i,:))-100;
    if ~isempty(s1)
        % decide the colour of the dots
        if touch_idx(i)==1
            plot(s1,i*ones(size(s1,2)),'.','Color',[1 0 0]),
        elseif touch_idx(i)==2
            plot(s1,i*ones(size(s1,2)),'.','Color',[0.5 0 0]),
        elseif touch_idx(i)>2
            plot(s1,i*ones(size(s1,2)),'.','Color',[0 0 0]),
        end
    end
end

plot([0 0],[0 size(spikes,1)],'Color',[0.5 0.5 0.5])
box off
xlabel('Time from touch onset [ms]')
ylabel('Touch number')
xlim([-100 100])
ylim([0 size(spikes,1)])
title('Order within trials')
hold off



filteredFR1=conv(mean(spikes(touch_idx==1,:)),boxcar,'same');
filteredFR2=conv(mean(spikes(touch_idx==2,:)),boxcar,'same');
filteredFR3=conv(mean(spikes(touch_idx==3,:)),boxcar,'same');

subplot(4,3,nplot(4))
hold on
plot(xtime,filteredFR1,'Color',[1 0 0])
plot(xtime,filteredFR2,'Color',[0.5 0 0])
plot(xtime,filteredFR3,'Color',[0 0 0])
box off
legend('First touch','Second touch','Later')
xlabel('Time from touch onset [ms]')
xlim([-100 100])
ylim([0 2])
plot([0 0],[0 2],'Color',[0.5 0.5 0.5])
hold off
end