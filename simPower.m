function StatPower=simPower(data4stats,idxIndependantEffect, rmModalities, indEffectNames, rmEffectNames)

% Define population parameters
numGroups = size(data4stats,2);
populationMeans = mean(data4stats);
populationSD = std(data4stats);
sampleSize = size(data4stats,1);
numSimulations = 1000;

% Initialize variables to store results
significantResults = zeros(1, numSimulations);

for i = 1:numSimulations
    % Generate simulated data
    data =[];
    for group = 1:numGroups-1
        data(:,group) = normrnd(populationMeans(group), populationSD(group), [1, sampleSize]);
    end

    % Perform ANOVA
    pValue = anova1(data, [], 'off');  % 'off' suppresses display

    % Check for significance (you may need to adjust the significance level)
    significance = pValue < 0.05;

    % Store result
    significantResults(i) = significance;
end

% Calculate power
power = mean(significantResults);



end