function [ES, diff, relDiff]=esCalculation0D(DATA)

n1 = numel(DATA{1});
n2 = numel(DATA{2});
nu = n1+n2-2;
Gnu = gamma(nu/2)/(sqrt(nu/2)*gamma((nu-1)/2));
SS1 = nanstd(DATA{1}).^2*(n1-1);
SS2 = nanstd(DATA{2}).^2*(n2-1);
pooledsd = sqrt((SS1 + SS2)/nu);
d = (nanmean(DATA{1})-nanmean(DATA{2}))./pooledsd;
ES =  abs(d);

diff=nanmean(DATA{1})-nanmean(DATA{2});
relDiff=100*(nanmean(DATA{1})-nanmean(DATA{2}))/nanmean(DATA{2});

end