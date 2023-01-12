function phFinal=formatPH4(ph,modalities,condNames)

string2compare=modalities{2};
for i=1:numel(string2compare)
    for j=1:size(ph,1)

        if all([contains(ph{j,2},string2compare{i}) contains(ph{j,3},string2compare{i})])
            keepRow(j,i)=i;
        else
            keepRow(j,i)=0;
        end

    end
end

keepRow=max(keepRow,[],2);
kRow=keepRow(keepRow>0);
ph2=ph(keepRow>0,:);

for k=1:2
    string2compare=modalities{2+k};
    for i=1:numel(string2compare)
        for j=1:size(ph2,1)

            if contains(ph2{j,1},string2compare{i})
                phT{j,k}=string2compare{i};
            end

        end
    end
end

string2compare=modalities{2};
for i=1:numel(string2compare)
    for j=1:size(ph2,1)

        if contains(ph2{j,2},string2compare{i})
            phT{j,3}=string2compare{i};
        end

    end
end

string2compare=modalities{1};
for i=1:numel(string2compare)
    for j=1:size(ph2,1)

        for k=1:2
            if contains(ph2{j,k+1},string2compare{i})
                phT{j,k+3}=string2compare{i};
            end

        end
    end
end

phT=[table(phT(:,1),'variableNames',condNames(3))  table(phT(:,2),'variableNames',condNames(4)) table(phT(:,3),'variableNames',condNames(2)) table(phT(:,4),'variableNames',{[condNames{1} '_1']}) table(phT(:,5),'variableNames',{[condNames{1} '_2']})];
phFinal=[phT, ph2(:,4:end)];

end