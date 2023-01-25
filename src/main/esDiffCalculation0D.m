function [ES, ESsd, diff, diffSD, relDiff, relDiffSD]=esDiffCalculation0D(DATA,isInd)

n1 = numel(DATA{1});
n2 = numel(DATA{2});
nu = n1+n2-2;
SS1 = nanstd(DATA{1}).^2*(n1-1);
SS2 = nanstd(DATA{2}).^2*(n2-1);
pooledsd = sqrt((SS1 + SS2)/nu);

if isInd
diff=nanmean(DATA{1})-nanmean(DATA{2});
relDiff=100*diff/nanmean(DATA{2});
diffSD=nan;
relDiffSD=nan;
else
diff=nanmean(DATA{1}-DATA{2});
relDiff=100*diff/nanmean(DATA{2});
diffSD=nanstd(DATA{1}-DATA{2});
relDiffSD=100*diffSD/nanmean(DATA{2});
end

ES = abs(diff)/pooledsd;
ESsd=sqrt((n1+n2)/(n1*n2)+(ES.^2)/(2*(n1+n2)));


end