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