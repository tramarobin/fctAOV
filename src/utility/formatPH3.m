function phFinal=formatPH3(ph,modalities,condNames)

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

ph3=[];
for k=1:max(double(ph2{:,1}))
    cat1=double(ph2{:,1});
    cut1=cat1==k;
    pht=ph2(cut1,:);
    kRt=kRow(cut1);

    for i=1:size(pht,1)
        n2{i,1}=modalities{2}{kRt(i)};
    end
    [a,b]=sort(n2);
    pht=pht(b,:);

    [whichPH]=findPHint(modalities{1},modalities{2});
    ph3=[ph3; pht(whichPH,:)];
end

clear keepRow
for i=1:numel(string2compare)
    for j=1:size(ph3,1)
        if all([contains(ph3{j,2},string2compare{i}) contains(ph3{j,3},string2compare{i})])
            keepRow(j,i)=i;
        else
            keepRow(j,i)=0;
        end
    end
end
keepRow=max(keepRow,[],2);
kRow=keepRow(keepRow>0);

for i=1:size(ph3,1)
    E2{i,1}=modalities{2}{kRow(i)};
    E1{i,1}=modalities{3}{ph3{i,1}};
    for e=1:2
        E31{i,e}=strrep(ph3{i,e+1},E2{i,1},'');
    end
end

varNames=string({condNames{3},condNames{2},[condNames{1} '_1'], [condNames{1} '_2']});
phT=table(E1, E2, E31(:,1), E31(:,2),VariableNames=varNames);
phFinal=[phT, ph3(:,4:end)];

end