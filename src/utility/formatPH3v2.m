function phFinal=formatPH3v2(ph,modalities,condNames)

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
nComp=0;
for i=1:numel(modalities{2})
    for j=1:numel(modalities{2})
        if i~=j
            nComp=nComp+1;
        end
    end
end
cutPs=1:nComp:numel(kRow);
cutPe=nComp:nComp:numel(kRow);
for k=1:numel(cutPs)
    ph2(cutPs(k):cutPe(k),:)=ph2(cutPs(k)-1+sort(kRow(cutPs(k):cutPe(k))),:);
end


string2compare=modalities{2};
for i=1:numel(string2compare)
    for j=1:size(ph2,1)

        if contains(ph2{j,2},string2compare{i})
            phT{j,1}=string2compare{i};
        end

    end
end

string2compare=modalities{1};
for i=1:numel(string2compare)
    for j=1:size(ph2,1)

        for k=1:2
            if contains(ph2{j,k+1},string2compare{i})
                phT{j,k+1}=string2compare{i};
            end

        end
    end
end


phT=[table(phT(:,1),'variableNames',condNames(2)) table(phT(:,2),'variableNames',{[condNames{1} '_1']}) table(phT(:,3),'variableNames',{[condNames{1} '_2']})];
phFinal=[ph2(:,1), phT, ph2(:,4:end)];

end