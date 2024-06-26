function Tables=AOVreformTable(Tables,stats, allMod, aov)

AVGt=[];
rmEffects=[{'allData'}, stats.condNames(find(stats.isRM))];
indEffects=stats.condNames(find(~stats.isRM));
nInd=numel(indEffects);
nRM=numel(rmEffects);

if ~isfield(stats,'effectOrder')
    stats.effectOrder=1:numel(stats.isRM);
end

for n=1:nRM

    effectRM=verifFieldName(rmEffects{n});

    dataEffect=Tables.(effectRM){1:numel(stats.ID),nInd+2:end};

    if ~isempty(findcol(fieldnames(stats),"addUnivariate"))
        dataEffect(:,end)=[];
    end

    AVG=mean(dataEffect,2,"omitnan");
    AVGall=mean(AVG,'all',"omitnan");
    SDall=std(AVG,[],"all","omitnan");

    AVGt=[AVG; nan; AVGall; SDall; nan];

    for i=1:nInd

        mods=allMod.(indEffects{i});

        for nMod=1:numel(mods)

            tMods=Tables.(effectRM){1:numel(stats.ID),1+i};
            p=findcol(tMods, mods{nMod});

            AVGp=mean(AVG(p));
            SDp=std(AVG(p),[],"all","omitnan");

            AVGt=[AVGt; AVGp; SDp; nan];

        end

    end


    Tables.(effectRM){:,size(Tables.(effectRM),2)+1}=AVGt;
    Tables.(effectRM).Properties.VariableNames{end}='ALL RM';





end


end
