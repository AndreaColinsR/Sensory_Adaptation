function balanced_anova(correct_1,incorrect_1,correct_later,incorrect_later,nsample)
%% balanced_anova resamples the data so that all entries have the same
% number of samples (nsample) 

nrep=50;

correct=[ones(nsample,1);zeros(nsample,1);ones(nsample,1);zeros(nsample,1)];
first_later=[ones(nsample,1);ones(nsample,1);ones(nsample,1)*2;ones(nsample,1)*2];

p=nan(2,nrep);

for i=1:nrep

y=[correct_1(randperm(numel(correct_1),nsample));...
    incorrect_1(randperm(numel(incorrect_1),nsample));...
    correct_later((randperm(numel(correct_later),nsample)));...
    incorrect_later((randperm(numel(incorrect_later),nsample)))];

p(:,i) = anovan(y,{correct,first_later},'display','off');
clear y
end

disp('-----------------------------------')
disp(['Balanced anova (' num2str(nrep) ' resamplings), N samples = ' num2str(nsample)])
disp (['p-value hit vs miss = ' num2str(mean(p(1,:))) ' +- ' num2str(std(p(1,:)))])
disp (['p-value first vs later = ' num2str(mean(p(2,:))) ' +- ' num2str(std(p(2,:)))])
end