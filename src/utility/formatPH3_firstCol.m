function phFinal=formatPH3firstCol(ph,modalities,condNames)

for k=1:2
    string2compare=modalities{k};
    for i=1:numel(string2compare)
        for j=1:size(ph,1)

            if contains(ph{j,1},string2compare{i})
                phT{j,k}=string2compare{i};
            end

        end
    end
end

phT=[table(phT(:,1),'variableNames',condNames(1)) table(phT(:,2),'variableNames',condNames(2))];
phFinal=[phT, ph(:,2:end)];

end