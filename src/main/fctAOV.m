% create graph and tables for a 4way anova with 2independant effect and 2rm depending on the
% project parameters values

function [Tables]=fctAOV(data,stats,display,units,saveDir)

%% Unpack structures and create save directory
warning('off');
if ~isempty(saveDir)
    mkdir(saveDir)
end
unpackStruct(stats)
unpackStruct(display)

%% Statistical names and order based on the names
% modalities
for nCond=1:numel(condNames)
    condNamesVerif{nCond}=verifFieldName(condNames{nCond});
end
loopRm=1;
for nCond=1:numel(cond4effect)
    allModalities{nCond}=unique(cond4effect{nCond},'stable');
    allModalitiesVerif{nCond}=verifFieldName(allModalities{nCond});
    allMod.(condNamesVerif{nCond})=unique(cond4effect{nCond},'stable');
end

% independant effect
indEffect=find(isRM==0);
idxIndependantEffect=zeros(size(data,1),numel(indEffect));
for nInd=1:numel(indEffect)
    modalitiesInd{nInd}=unique(cond4effect{indEffect(nInd)},'stable');
    for m=1:numel(modalitiesInd{nInd})
        idxIndEffect=findcol(cond4effect{indEffect(nInd)},modalitiesInd{nInd}{m});
        idxIndependantEffect(idxIndEffect,nInd)=m;
    end
end

% rm effect
rmEffect=find(isRM==1);
for nRm=1:sum(isRM)
    for nModRm=1:numel(cond4effect{rmEffect(nRm)})
        modalitiesRM{nRm}=unique(cond4effect{rmEffect(nRm)},'stable');
    end
end

if sum(isRM)==0
    condNames4Table={'No RM'};
elseif sum(isRM)==1
    effectRMaov{1}=transpose(string(cond4effect{rmEffect(1)}));
    effectRM{1}=cond4effect{rmEffect(1)};
    condNames4Table=cond4effect{rmEffect(1)};
elseif sum(isRM)==2
    effectRMaov{1}=transpose(repmat(string(cond4effect{rmEffect(1)}),1,numel(cond4effect{rmEffect(2)})));
    for i=1:numel(cond4effect{rmEffect(1)})
        for j=1:numel(cond4effect{rmEffect(2)})
            effectRM{1}{numel(cond4effect{rmEffect(2)})*(i-1)+j}=cond4effect{rmEffect(1)}{i};
            condNames4Table{numel(cond4effect{rmEffect(2)})*(i-1)+j}=[cond4effect{rmEffect(1)}{i} ' ' cond4effect{rmEffect(2)}{j}];
        end
    end
    effectRMaov{2}=[];
    for j=1:numel(cond4effect{rmEffect(2)})
        for i=1:numel(cond4effect{rmEffect(1)})
            effectRMaov{2}=[effectRMaov{2}; string(cond4effect{rmEffect(2)}{j})];
        end
    end
    effectRM{2}=[repmat(cond4effect{rmEffect(2)},1,numel(cond4effect{rmEffect(1)}))];
end

% order 2 selct means for ES calculations
if sum(isRM)==1
    col4means{1}=transpose(1:size(data,2));
elseif sum(isRM)==2
    for i=1:numel(cond4effect{rmEffect(1)})
        col4means{1}(i,:)=numel(cond4effect{rmEffect(2)})*(i-1)+1:numel(cond4effect{rmEffect(2)})*(i-1)+numel(cond4effect{rmEffect(2)});
    end
    for j=1:numel(cond4effect{rmEffect(2)})
        col4means{2}(j,:)=j:numel(cond4effect{rmEffect(2)}):numel(condNames4Table);
    end
end

if sum(isRM)>1
    for nCond1=1:sum(isRM)
        for nCond2=1:sum(isRM)
            if nCond1~=nCond2
                order4EStmp=[];
                for k=1:size(col4means{nCond1},2)
                    loop=1;
                    comp=col4means{nCond1}(:,k);
                    for i=1:numel(comp)
                        for j=1:numel(comp)
                            if comp(i)<comp(j)
                                order(loop,1)=comp(i);
                                order(loop,2)=comp(j);
                                loop=loop+1;
                            end
                        end
                    end
                    [order(:,1),b]=sort(order(:,1),'descend');
                    order(:,2)=order(b,2);
                    order4EStmp=[order4EStmp; order];
                    clear order
                end
                order4ES.([condNamesVerif{rmEffect(nCond1)} 'By' condNamesVerif{rmEffect(nCond2)}])=order4EStmp;
            end
        end
    end
end

% Change the configuration of data
numSub=size(data,1);
if sum(isRM)==1
    if iscell(data)
        for i=1:numel(cond4effect{rmEffect(1)})
            for s=1:size(data,1)
                varData(s,i)=nanmean(data{s,i});
            end
            t{i}=varData(:,i);
        end
    else
        varData=data;
        for i=1:numel(cond4effect{rmEffect(1)})
            t{i}=varData(:,i);
        end
    end
    clear varData
elseif sum(isRM)==2
    for i=1:numel(cond4effect{rmEffect(1)})
        for s=1:size(data,1)
            for j=1:numel(cond4effect{rmEffect(2)})
                if iscell(data)
                    varData(s,j)=nanmean(data{s,i,j});
                else
                    varData(s,j)=data(s,i,j);
                end
            end
        end
        t{i}=varData;
        clear varData
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create tables for means
if numel(indEffect)==0
    cNames=['ID' condNames4Table];
elseif numel(indEffect)==1
    cNames=['ID' condNamesVerif{indEffect(1)} condNames4Table];
elseif numel(indEffect)==2
    cNames=['ID' condNamesVerif{indEffect(1)} condNamesVerif{indEffect(2)} condNames4Table];
end

tXl=[];
if ~isempty(nRm)
    for i=1:numel(cond4effect{rmEffect(1)})
        tXl=[tXl t{i}];
    end
else
    if iscell(data)
        for s=1:size(data,1)
            tXl(s,1)=nanmean(data{s});
        end
    else
        tXl=data;
    end

end
% nan the subjects with issues
% tXl(isnan(mean(tXl,2)),:)=nan;
tXL(:,1)=table([ID(1:numSub); nan]);
for nInd=1:numel(indEffect)
    tXL(:,nInd+1)=table([cond4effect{indEffect(nInd)}(1:numSub); nan]);
end
tXLmeansCol={'All participants mean'; 'All participants sd';' '};
for nInd=1:numel(indEffect)
    for nMod=1:numel(modalitiesInd{nInd})
        tXLmeansCol=[tXLmeansCol; {[modalitiesInd{nInd}{nMod} ' mean']; [modalitiesInd{nInd}{nMod} ' sd']; ' '}];
    end
end
if numel(indEffect)>1
    for nInd1=1:numel(indEffect)
        for nInd2=1:numel(indEffect)
            if nInd2>nInd1
                for nMod1=1:numel(modalitiesInd{nInd1})
                    for nMod2=1:numel(modalitiesInd{nInd2})
                        tXLmeansCol=[tXLmeansCol; {[modalitiesInd{nInd1}{nMod1} ' & ' modalitiesInd{nInd2}{nMod2}  ' mean']; [modalitiesInd{nInd1}{nMod1} ' & ' modalitiesInd{nInd2}{nMod2}  ' sd']; ' '}];
                    end
                end
            end
        end
    end
end
if numel(indEffect)>2
    for nInd1=1:numel(indEffect)
        for nInd2=1:numel(indEffect)
            for nInd3=1:numel(indEffect)
                if nInd3>nInd2
                    if nInd2>nInd1
                        for nMod1=1:numel(modalitiesInd{nInd1})
                            for nMod2=1:numel(modalitiesInd{nInd2})
                                for nMod3=1:numel(modalitiesInd{nInd3})
                                    tXLmeansCol=[tXLmeansCol; {[modalitiesInd{nInd1}{nMod1} ' & ' modalitiesInd{nInd2}{nMod2} ' & ' modalitiesInd{nInd3}{nMod3}  ' mean']; [modalitiesInd{nInd1}{nMod1} ' & ' modalitiesInd{nInd2}{nMod2} ' & ' modalitiesInd{nInd3}{nMod3}  ' sd']; ' '}];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

tXLmeans=table(tXLmeansCol);
tXLmeans(:,2:numel(indEffect)+1)=table({nan});

for i=1:numel(condNames4Table)
    tXL(:,numel(indEffect)+1+i)=table([tXl(:,i); nan]);
    tXLmeans(1,numel(indEffect)+1+i)=table(nanmean(tXl(:,i)));
    tXLmeans(2,numel(indEffect)+1+i)=table(nanstd(tXl(:,i)));
    tXLmeans(3:end,numel(indEffect)+1+i)=table(nan);
end
data4plot.allData=tXl;
if ~isempty(findcol(fieldnames(stats),"addUnivariate"))
    data4plot.allData(:,end)=[];
end
loop=4;
if numel(indEffect)>0
    for nInd1=1:numel(indEffect)
        for nMod1=1:numel(modalitiesInd{nInd1})
            for i=1:numel(condNames4Table)
                tXLmeans(loop,numel(indEffect)+1+i)=table(nanmean(tXl(idxIndependantEffect(:,nInd1)==nMod1,i)));
                tXLmeans(loop+1,numel(indEffect)+1+i)=table(nanstd(tXl(idxIndependantEffect(:,nInd1)==nMod1,i)));
                tXLmeans(loop+2,numel(indEffect)+1+i)=table(nan);
            end
            data4plot.(condNamesVerif{[indEffect(nInd1)]}).(allModalitiesVerif{[indEffect(nInd1)]}{nMod1})=data4plot.allData(idxIndependantEffect(:,nInd1)==nMod1,:);
            loop=loop+3;
        end
    end
end
if numel(indEffect)>1
    for nInd1=1:numel(indEffect)
        for nInd2=1:numel(indEffect)
            if nInd2>nInd1
                for nMod1=1:numel(modalitiesInd{nInd1})
                    for nMod2=1:numel(modalitiesInd{nInd2})
                        for i=1:numel(condNames4Table)
                            tXLmeans(loop,numel(indEffect)+1+i)=table(nanmean(tXl(idxIndependantEffect(:,nInd1)==nMod1 & idxIndependantEffect(:,nInd2)==nMod2,i)));
                            tXLmeans(loop+1,numel(indEffect)+1+i)=table(nanstd(tXl(idxIndependantEffect(:,nInd1)==nMod1 & idxIndependantEffect(:,nInd2)==nMod2,i)));
                            tXLmeans(loop+2,numel(indEffect)+1+i)=table(nan);
                        end
                        data4plot.([condNamesVerif{[indEffect(nInd1)]} 'By' condNamesVerif{[indEffect(nInd2)]}]).([allModalitiesVerif{[indEffect(nInd1)]}{nMod1} allModalitiesVerif{[indEffect(nInd2)]}{nMod2}])=data4plot.allData(idxIndependantEffect(:,nInd1)==nMod1 & idxIndependantEffect(:,nInd2)==nMod2,:);
                        data4plot.([condNamesVerif{[indEffect(nInd1)]} 'By' condNamesVerif{[indEffect(nInd2)]}]).(allModalitiesVerif{[indEffect(nInd1)]}{nMod1}).(allModalitiesVerif{indEffect(nInd2)}{nMod2})=data4plot.allData(idxIndependantEffect(:,nInd1)==nMod1 & idxIndependantEffect(:,nInd2)==nMod2,:);
                        data4plot.([condNamesVerif{[indEffect(nInd2)]} 'By' condNamesVerif{[indEffect(nInd1)]}]).(allModalitiesVerif{[indEffect(nInd2)]}{nMod2}).(allModalitiesVerif{indEffect(nInd1)}{nMod1})=data4plot.allData(idxIndependantEffect(:,nInd1)==nMod1 & idxIndependantEffect(:,nInd2)==nMod2,:);
                        loop=loop+3;
                    end
                end
            end
        end
    end
end
if numel(indEffect)>2
    for nInd1=1:numel(indEffect)
        for nInd2=1:numel(indEffect)
            for nInd3=1:numel(indEffect)
                if nInd3>nInd2
                    if nInd2>nInd1
                        for nMod1=1:numel(modalitiesInd{nInd1})
                            for nMod2=1:numel(modalitiesInd{nInd2})
                                for nMod3=1:numel(modalitiesInd{nInd3})
                                    for i=1:numel(condNames4Table)
                                        tXLmeans(loop,numel(indEffect)+1+i)=table(nanmean(tXl(idxIndependantEffect(:,nInd1)==nMod1 & idxIndependantEffect(:,nInd2)==nMod2 & idxIndependantEffect(:,nInd3)==nMod3,i)));
                                        tXLmeans(loop+1,numel(indEffect)+1+i)=table(nanstd(tXl(idxIndependantEffect(:,nInd1)==nMod1 & idxIndependantEffect(:,nInd2)==nMod2 & idxIndependantEffect(:,nInd3)==nMod3,i)));
                                        tXLmeans(loop+2,numel(indEffect)+1+i)=table(nan);
                                    end
                                    loop=loop+3;
                                    data4plot.([condNamesVerif{[indEffect(nInd1)]} 'By' condNamesVerif{[indEffect(nInd2)]} 'By' condNamesVerif{[indEffect(nInd3)]}]).([allModalitiesVerif{[indEffect(nInd1)]}{nMod1} allModalitiesVerif{[indEffect(nInd2)]}{nMod2} allModalitiesVerif{[indEffect(nInd3)]}{nMod3}])=data4plot.allData(idxIndependantEffect(:,nInd1)==nMod1 & idxIndependantEffect(:,nInd2)==nMod2 & idxIndependantEffect(:,nInd3)==nMod3,:);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

tXL.Properties.VariableNames=cNames;
tXLmeans.Properties.VariableNames=cNames;
tablemeans.allData=[tXL; tXLmeans];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Stats
% aov and posthoc
if ~isempty(nRm)
    for k=1:numel(cond4effect{rmEffect(1)})
        data4stats(:,k,:)=t{k};
    end
else
    data4stats(:,1,:)=tXl;
end

if numel(rmEffect)==0 & numel(indEffect)==1
    [p,tbl,statsAOV]=anova1(data4stats,idxIndependantEffect,'off');
end

if numel(rmEffect)==1 & numel(indEffect)==0
    [tbl,rm]=mixed_anova(data4stats,[], {condNamesVerif{rmEffect(1)}});
    rm.WithinDesign.(condNamesVerif{rmEffect(1)})=effectRMaov{1};
end

if numel(rmEffect)==2 & numel(indEffect)==0
    [tbl,rm]=mixed_anova(data4stats,[], {condNamesVerif{rmEffect(1)} condNamesVerif{rmEffect(2)}});
    rm.WithinDesign.(condNamesVerif{rmEffect(1)})=effectRMaov{1};
    rm.WithinDesign.(condNamesVerif{rmEffect(2)})=effectRMaov{2};
end

if numel(rmEffect)==2 & numel(indEffect)==2

    [tbl,rm]=mixed_anova(data4stats,idxIndependantEffect, {condNamesVerif{rmEffect(1)}, condNamesVerif{rmEffect(2)}}, {condNamesVerif{indEffect(1)}, condNamesVerif{indEffect(2)}});
    rm.WithinDesign.(condNamesVerif{rmEffect(1)})=effectRMaov{1};
    rm.WithinDesign.(condNamesVerif{rmEffect(2)})=effectRMaov{2};

elseif numel(rmEffect)==1 & numel(indEffect)==1

    [tbl,rm]=mixed_anova(data4stats,idxIndependantEffect, {condNamesVerif{rmEffect(1)}}, {condNamesVerif{indEffect(1)}});
    rm.WithinDesign.(condNamesVerif{rmEffect(1)})=effectRMaov{1};

elseif numel(rmEffect)==1 & numel(indEffect)==2

    [tbl,rm]=mixed_anova(data4stats,idxIndependantEffect, {condNamesVerif{rmEffect(1)}}, {condNamesVerif{indEffect(1)}, condNamesVerif{indEffect(2)}});
    rm.WithinDesign.(condNamesVerif{rmEffect(1)})=effectRMaov{1};


elseif numel(rmEffect)==2 & numel(indEffect)==1

    [tbl,rm]=mixed_anova(data4stats,idxIndependantEffect, {condNamesVerif{rmEffect(1)}, condNamesVerif{rmEffect(2)}}, {condNamesVerif{indEffect(1)}});
    rm.WithinDesign.(condNamesVerif{rmEffect(1)})=effectRMaov{1};
    rm.WithinDesign.(condNamesVerif{rmEffect(2)})=effectRMaov{2};

end

if istable(tbl)
    rNames=tbl.Properties.RowNames;
    cNamest=tbl.Properties.VariableNames;

    rCol=ones(size(tbl,1),1);
    rCol(1)=0;
    rCol(findcol(rNames,'Error'))=0;
    rCol=logical(rCol);
    rNames=rNames(rCol);
    for i=1:numel(rNames)
        rNames{i}=strrep(rNames{i}, '(Intercept):','');
    end
    for i=1:size(tbl,2)
        aovt(:,i)=table(tbl{rCol,i});
    end

    aovt.Properties.VariableNames=cNamest;
    aov=[table(rNames,'VariableNames',{'Effect'}) aovt];

else

    rNames=tbl(2:end,1);
    cNamest=tbl(1,:);

    aov=splitvars(table(tbl(2,:)));
    aov{1,1}=condNames(1);
    aov.Properties.VariableNames={'Effect', 'SumSq', 'DF', 'MeanSq', 'F', 'pValue'};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1 effect
if numel(rmEffect)==0 & numel(indEffect)==1

    ph=multcompare(statsAOV,"CriticalValueType",postHocType,"Display","off");
    [whichPH,order4ES.(condNamesVerif{1})]=findPH(allModalities{1});
    %     phNew=ph(whichPH,:);

    phNew(:,1)=table(modalitiesInd{1}(ph(:,1)));
    phNew(:,2)=table(modalitiesInd{1}(ph(:,2)));
    phNew(:,3)=table(ph(:,6));
    phNew.Properties.VariableNames={[condNamesVerif{1} '_1'], [condNamesVerif{1} '_2'], 'pValue'};
    postHoc.(condNamesVerif{1})=phNew;

else

    for nCond=1:numel(condNames)
        ph=multcompare(rm,condNamesVerif{nCond},'ComparisonType',postHocType);
        [whichPH,order4ES.(condNamesVerif{nCond})]=findPH(allModalities{nCond});
        postHoc.(condNamesVerif{nCond})=ph(whichPH,:);
        if any(nCond==indEffect)
            postHoc.(condNamesVerif{nCond}).([condNamesVerif{nCond} '_1'])=modalitiesInd{nCond}(postHoc.(condNamesVerif{nCond}).([condNamesVerif{nCond} '_1']));
            postHoc.(condNamesVerif{nCond}).([condNamesVerif{nCond} '_2'])=modalitiesInd{nCond}(postHoc.(condNamesVerif{nCond}).([condNamesVerif{nCond} '_2']));
        end
    end
end
%% 2 effects
if numel(condNames)>1
    for nCond1=1:numel(condNames)
        for nCond2=1:numel(condNames)
            if nCond2~=nCond1
                ph=multcompare(rm,condNamesVerif{nCond1},'By',condNamesVerif{nCond2},'ComparisonType',postHocType);
                [whichPH]=findPHint(allModalities{nCond1},allModalities{nCond2});
                postHoc.([condNamesVerif{nCond1} 'By'  condNamesVerif{nCond2}])=ph(whichPH,:);
                if any(nCond1==indEffect)
                    postHoc.([condNamesVerif{nCond1} 'By'  condNamesVerif{nCond2}]).([condNamesVerif{nCond1} '_1'])=modalitiesInd{nCond1}(postHoc.([condNamesVerif{nCond1} 'By'  condNamesVerif{nCond2}]).([condNamesVerif{nCond1} '_1']));
                    postHoc.([condNamesVerif{nCond1} 'By'  condNamesVerif{nCond2}]).([condNamesVerif{nCond1} '_2'])=modalitiesInd{nCond1}(postHoc.([condNamesVerif{nCond1} 'By'  condNamesVerif{nCond2}]).([condNamesVerif{nCond1} '_2']));
                end
                if any(nCond2==indEffect)
                    postHoc.([condNamesVerif{nCond1} 'By'  condNamesVerif{nCond2}]).([condNamesVerif{nCond2}])=modalitiesInd{nCond2}(postHoc.([condNamesVerif{nCond1} 'By'  condNamesVerif{nCond2}]).([condNamesVerif{nCond2}]));
                end
            end
        end
    end
end


%% 3 effects
if numel(condNames)>2
    for nCond1=1:numel(condNames)
        for nCond2=1:numel(condNames)
            for nCond3=1:numel(condNames)

                % ind/rm*rm*rm
                if all([nCond3~=nCond2 & nCond3~=nCond1 & nCond2~=nCond1 any(nCond2==rmEffect) any(nCond3==rmEffect)])

                    if numel(condNames)>3
                        [tbl,rm]=mixed_anova(tXl,idxIndependantEffect, {[condNamesVerif{nCond2} 'By' condNamesVerif{nCond3}]}, {condNamesVerif{indEffect(1)}, condNamesVerif{indEffect(2)}});
                    elseif numel(condNames)==3 & any(nCond1==indEffect)
                        [tbl,rm]=mixed_anova(tXl,idxIndependantEffect, {[condNamesVerif{nCond2} 'By' condNamesVerif{nCond3}]}, {condNamesVerif{indEffect(nCond1)}});
                    end

                    rm.WithinDesign.([condNamesVerif{nCond2} 'By' condNamesVerif{nCond3}])=transpose(string(condNames4Table));
                    ph=multcompare(rm,condNamesVerif{nCond1},'By', [condNamesVerif{nCond2} 'By' condNamesVerif{nCond3}], 'ComparisonType',postHocType);
                    [whichPH]=findPHint(allModalities{nCond1},condNames4Table);
                    if any(nCond1==indEffect)
                        ph.([condNamesVerif{nCond1} '_1'])=modalitiesInd{nCond1}(ph.([condNamesVerif{nCond1} '_1']));
                        ph.([condNamesVerif{nCond1} '_2'])=modalitiesInd{nCond1}(ph.([condNamesVerif{nCond1} '_2']));
                    end
                    ph=ph(whichPH,:);
                    phNew=formatPH3_firstCol(ph,allModalities([nCond2 nCond3]), condNames([nCond2 nCond3]));
                    postHoc.([condNamesVerif{nCond1} 'By'  condNamesVerif{nCond2} 'By'  condNamesVerif{nCond3}])=phNew;

                    % ind*ind*ind/rm
                elseif all([nCond3~=nCond2 & nCond3~=nCond1 & nCond2~=nCond1 any(nCond1==indEffect) any(nCond2==indEffect)])

                    for s=1:size(idxIndependantEffect,1)
                        int2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}])(s,1)=string([cond4effect{nCond1}{s} '' cond4effect{nCond2}{s}]);
                    end
                    modInt2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}])=unique(int2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}]));
                    allMod.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}])=unique(cellstr(int2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}])));
                    for s=1:numel(int2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}]))
                        indInt2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}])(s,1)=findcol(modInt2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}]), int2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}])(s));
                    end

                    [tbl,rm]=mixed_anova(data4stats,indInt2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}]), {condNamesVerif{nCond3}}, {[condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}]});
                    if any(nCond3==rmEffect)
                        rm.WithinDesign.(condNamesVerif{rmEffect(find(nCond3==rmEffect))})=effectRMaov{find(nCond3==rmEffect)};
                    end

                    ph=multcompare(rm,[condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}], 'By', condNamesVerif{nCond3}, 'ComparisonType',postHocType);
                    [whichPH]=findPHint(modInt2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}]),allModalities{nCond3});
                    [~, order4ES.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}])]=findPH(modInt2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}]));
                    ph=ph(whichPH,:);
                    ph.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2} '_1'])=modInt2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}])(ph.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2} '_1']));
                    ph.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2} '_2'])=modInt2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}])(ph.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2} '_2']));
                    phNew=formatPH3v2(ph,allModalities([nCond1 nCond2 nCond3]), condNames([nCond1 nCond2 nCond3]));
                    postHoc.([condNamesVerif{nCond1} 'By'  condNamesVerif{nCond2} 'By'  condNamesVerif{nCond3}])=phNew;

                    %ind/rm*ind*ind
                elseif all([nCond3>nCond2 & nCond3~=nCond1 & nCond2~=nCond1 any(nCond2==indEffect) any(nCond3==indEffect)])

                    for s=1:size(idxIndependantEffect,1)
                        int2.([condNamesVerif{nCond2} 'By' condNamesVerif{nCond3}])(s,1)=string([cond4effect{nCond2}{s} '' cond4effect{nCond3}{s}]);
                    end
                    modInt2.([condNamesVerif{nCond2} 'By' condNamesVerif{nCond3}])=unique(int2.([condNamesVerif{nCond2} 'By' condNamesVerif{nCond3}]));
                    for s=1:numel( int2.([condNamesVerif{nCond2} 'By' condNamesVerif{nCond3}]))
                        indInt2.([condNamesVerif{nCond2} 'By' condNamesVerif{nCond3}])(s,1)=findcol(modInt2.([condNamesVerif{nCond2} 'By' condNamesVerif{nCond3}]), int2.([condNamesVerif{nCond2} 'By' condNamesVerif{nCond3}])(s));
                    end
                    if any(nCond1==indEffect)
                        if numel(condNames)>3
                            [tbl,rm]=mixed_anova(data4stats,indInt2.([condNamesVerif{nCond2} 'By' condNamesVerif{nCond3}]), {condNamesVerif{rmEffect(1)}, condNamesVerif{rmEffect(2)}}, {[condNamesVerif{nCond2} 'By' condNamesVerif{nCond3}]});
                            rm.WithinDesign.(condNamesVerif{rmEffect(1)})=effectRMaov{1};
                            rm.WithinDesign.(condNamesVerif{rmEffect(2)})=effectRMaov{2};
                        end
                    elseif any(nCond1==rmEffect)
                        if numel(condNames)==3
                            [tbl,rm]=mixed_anova(data4stats,indInt2.([condNamesVerif{nCond2} 'By' condNamesVerif{nCond3}]), {condNamesVerif{nCond1}}, {[condNamesVerif{nCond2} 'By' condNamesVerif{nCond3}]});
                            rm.WithinDesign.(condNamesVerif{nCond1})=effectRMaov{find(nCond1==rmEffect)};

                        else
                            [tbl,rm]=mixed_anova(data4stats,indInt2.([condNamesVerif{nCond2} 'By' condNamesVerif{nCond3}]), {condNamesVerif{rmEffect(1)}, condNamesVerif{rmEffect(2)}}, {[condNamesVerif{nCond2} 'By' condNamesVerif{nCond3}]});
                            rm.WithinDesign.(condNamesVerif{rmEffect(1)})=effectRMaov{1};
                            rm.WithinDesign.(condNamesVerif{rmEffect(2)})=effectRMaov{2};
                        end
                    end
                    ph=multcompare(rm,condNamesVerif{nCond1},'By', [condNamesVerif{nCond2} 'By' condNamesVerif{nCond3}], 'ComparisonType',postHocType);
                    [whichPH]=findPHint(allModalities{nCond1},modInt2.([condNamesVerif{nCond2} 'By' condNamesVerif{nCond3}]));
                    ph=ph(whichPH,:);
                    ph.([condNamesVerif{nCond2} 'By' condNamesVerif{nCond3}])=modInt2.([condNamesVerif{nCond2} 'By' condNamesVerif{nCond3}])(ph.([condNamesVerif{nCond2} 'By' condNamesVerif{nCond3}]));
                    phNew=formatPH3_firstCol(ph,allModalities([nCond2 nCond3]), condNames([nCond2 nCond3]));
                    postHoc.([condNamesVerif{nCond1} 'By'  condNamesVerif{nCond2} 'By'  condNamesVerif{nCond3}])=phNew;

                    % rm*rm*ind
                elseif all([nCond3~=nCond2 & nCond3~=nCond1 & nCond2~=nCond1 any(nCond1==rmEffect) any(nCond2==rmEffect)]) & any(nCond3==indEffect)

                    [tbl,rm]=mixed_anova(tXl,idxIndependantEffect(:,indEffect(nCond3)), {[condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}]}, {condNamesVerif{indEffect(nCond3)}});
                    rm.WithinDesign.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}])=condNames4Table';
                    ph=multcompare(rm, [condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}], 'By', condNamesVerif{nCond3}, 'ComparisonType',postHocType);
                    phNew=formatPH3(ph,allModalities([nCond1 nCond2 nCond3]), condNames([nCond1 nCond2 nCond3]));
                    postHoc.([condNamesVerif{nCond1} 'By'  condNamesVerif{nCond2} 'By'  condNamesVerif{nCond3}])=phNew;

                end
            end
        end
    end
end

%% 4 effects
% ind effect first
if numel(condNames)>3
    for nCond1=1:numel(condNames)
        for nCond2=1:numel(condNames)
            for nCond3=1:numel(condNames)
                for nCond4=1:numel(condNames)

                    if all([nCond1~=nCond2 & nCond1~=nCond3 & nCond1~=nCond4 & nCond4>nCond3 & nCond2~=nCond3 & nCond2~=nCond4 any(nCond1==indEffect) any(nCond2==indEffect) any(nCond3==rmEffect) any(nCond4==rmEffect)])

                        for s=1:size(idxIndependantEffect,1)
                            int2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}])(s,1)=string([cond4effect{nCond1}{s} ' & ' cond4effect{nCond2}{s}]);
                        end
                        modInt2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}])=unique(int2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}]));
                        for s=1:numel(int2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}]))
                            indInt2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}])(s,1)=findcol(modInt2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}]),int2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}])(s));
                        end
                        [tbl,rm]=mixed_anova(tXl,indInt2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}]), {[condNamesVerif{nCond3} 'By' condNamesVerif{nCond4}]}, {[condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}]});
                        rm.WithinDesign.([condNamesVerif{nCond3} 'By' condNamesVerif{nCond4}])=transpose(string(condNames4Table));
                        ph=multcompare(rm,[condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}],'By', [condNamesVerif{nCond3} 'By' condNamesVerif{nCond4}], 'ComparisonType',postHocType);
                        [whichPH]=findPHint(modInt2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}]),condNames4Table);
                        phName=[condNamesVerif{nCond1}  'By'  condNamesVerif{nCond2} 'By' condNamesVerif{nCond3} 'By'  condNamesVerif{nCond4}];
                        ind2replaceName=[condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}];
                        ph=ph(whichPH,:);
                        ph.([ind2replaceName '_1'])=modInt2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}])(ph.([ind2replaceName '_1']));
                        ph.([ind2replaceName '_2'])=modInt2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}])(ph.([ind2replaceName '_2']));
                        phNew=formatPH4(ph,allModalities([nCond1 nCond2 nCond3 nCond4]), condNames([nCond1 nCond2 nCond3 nCond4]));
                        postHoc.(phName)=phNew;

                    end
                end
            end
        end
    end
end

% rm effect first
if numel(condNames)>3
    for nCond1=1:numel(condNames)
        for nCond2=1:numel(condNames)
            for nCond3=1:numel(condNames)
                for nCond4=1:numel(condNames)

                    if all([nCond2>nCond1 & nCond1~=nCond3 & nCond1~=nCond4 & nCond4~=nCond3 & nCond2~=nCond3 & nCond2~=nCond4 any(nCond1==indEffect) any(nCond2==indEffect) any(nCond3==rmEffect) any(nCond4==rmEffect)])

                        for s=1:size(idxIndependantEffect,1)
                            int2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}])(s,1)=string([cond4effect{nCond1}{s} ' & ' cond4effect{nCond2}{s}]);
                        end
                        modInt2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}])=unique(int2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}]));
                        for s=1:numel(int2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}]))
                            indInt2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}])(s,1)=findcol(modInt2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}]),int2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}])(s));
                        end

                        [tbl,rm]=mixed_anova(tXl,indInt2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}]), {[condNamesVerif{nCond3} 'By' condNamesVerif{nCond4}]}, {[condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}]});
                        rm.WithinDesign.([condNamesVerif{nCond3} 'By' condNamesVerif{nCond4}])=transpose(string(condNames4Table));

                        ph=multcompare(rm,[condNamesVerif{nCond3} 'By' condNamesVerif{nCond4}],'By', [condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}], 'ComparisonType',postHocType);
                        [whichPH]=findPHint(condNames4Table, modInt2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}]));
                        phName=[condNamesVerif{nCond3}  'By'  condNamesVerif{nCond4} 'By' condNamesVerif{nCond1} 'By'  condNamesVerif{nCond2}];
                        ind2replaceName=[condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}];
                        ph=ph(whichPH,:);
                        ph.(ind2replaceName)=modInt2.([condNamesVerif{nCond1} 'By' condNamesVerif{nCond2}])(ph.(ind2replaceName));

                        phNew=formatPH4(ph,allModalities([nCond3 nCond4 nCond1 nCond2]), condNames([nCond3 nCond4 nCond1 nCond2]));
                        postHoc.(phName)=phNew;

                    end
                end
            end
        end
    end
end

%% Tables for means + Effect sizes for rm effect
for nCond=1:numel(rmEffect)

    for b=1:size(col4means{nCond},1)
        means(:,b)=mean(tXl(:,col4means{nCond}(b,:)),2);
    end
    means(isnan(mean(means,2)),:)=nan;
    tempTable(:,1)=table([ID(1:numSub); nan]);
    for nInd=1:numel(indEffect)
        tempTable(:,nInd+1)=table([cond4effect{indEffect(nInd)}(1:numSub); nan]);
    end
    for i=1:numel(cond4effect{rmEffect(nCond)})
        tempTable(:,1+i+numel(indEffect))=table([means(:,i); nan]);
    end
    tempTableMeans=table(tXLmeansCol,VariableNames="Var1");
    tempTableMeans(:,2:numel(indEffect)+1)=table({nan});
    for i=1:numel(cond4effect{rmEffect(nCond)})
        tempTableMeans(1,numel(indEffect)+1+i)=table(nanmean(means(:,i)));
        tempTableMeans(2,numel(indEffect)+1+i)=table(nanstd(means(:,i)));
        tempTableMeans(3,numel(indEffect)+1+i)=table(nan);
    end

    % adding independant effect in tables
    loop=4;
    if numel(indEffect)>0
        for nInd1=1:numel(indEffect)
            for nMod1=1:numel(modalitiesInd{nInd1})
                for i=1:numel(cond4effect{rmEffect(nCond)})
                    tempTableMeans(loop,numel(indEffect)+1+i)=table(nanmean(means(idxIndependantEffect(:,nInd1)==nMod1,i)));
                    tempTableMeans(loop+1,numel(indEffect)+1+i)=table(nanstd(means(idxIndependantEffect(:,nInd1)==nMod1,i)));
                    tempTableMeans(loop+2,numel(indEffect)+1+i)=table(nan);
                end
                loop=loop+3;
            end
        end
    end
    if numel(indEffect)>1
        for nInd1=1:numel(indEffect)
            for nInd2=1:numel(indEffect)
                if nInd2>nInd1
                    for nMod1=1:numel(modalitiesInd{nInd1})
                        for nMod2=1:numel(modalitiesInd{nInd2})
                            for i=1:numel(cond4effect{rmEffect(nCond)})
                                tempTableMeans(loop,numel(indEffect)+1+i)=table(nanmean(means(idxIndependantEffect(:,nInd1)==nMod1 & idxIndependantEffect(:,nInd2)==nMod2,i)));
                                tempTableMeans(loop+1,numel(indEffect)+1+i)=table(nanstd(means(idxIndependantEffect(:,nInd1)==nMod1 & idxIndependantEffect(:,nInd2)==nMod2,i)));
                                tempTableMeans(loop+2,numel(indEffect)+1+i)=table(nan);
                            end
                            loop=loop+3;
                        end
                    end
                end
            end
        end
    end
    if numel(indEffect)>2
        for nInd1=1:numel(indEffect)
            for nInd2=1:numel(indEffect)
                for nInd3=1:numel(indEffect)
                    if nInd3>nInd2
                        if nInd2>nInd1
                            for nMod1=1:numel(modalitiesInd{nInd1})
                                for nMod2=1:numel(modalitiesInd{nInd2})
                                    for nMod3=1:numel(modalitiesInd{nInd3})
                                        for i=1:numel(cond4effect{rmEffect(nCond)})
                                            tempTableMeans(loop,numel(indEffect)+1+i)=table(nanmean(means(idxIndependantEffect(:,nInd1)==nMod1 & idxIndependantEffect(:,nInd2)==nMod2 & idxIndependantEffect(:,nInd3)==nMod3,i)));
                                            tempTableMeans(loop+1,numel(indEffect)+1+i)=table(nanstd(means(idxIndependantEffect(:,nInd1)==nMod1 & idxIndependantEffect(:,nInd2)==nMod2 & idxIndependantEffect(:,nInd3)==nMod3,i)));
                                            tempTableMeans(loop+2,numel(indEffect)+1+i)=table(nan);
                                        end
                                        loop=loop+3;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    tablemeans.(condNamesVerif{rmEffect(nCond)})=[tempTable; tempTableMeans];
    if numel(indEffect)==0
        tablemeans.(condNamesVerif{rmEffect(nCond)}).Properties.VariableNames=['ID' cond4effect{rmEffect(nCond)}];
    elseif numel(indEffect)==1
        tablemeans.(condNamesVerif{rmEffect(nCond)}).Properties.VariableNames=['ID' condNamesVerif{indEffect(1)} cond4effect{rmEffect(nCond)}];
    elseif numel(indEffect)==2
        tablemeans.(condNamesVerif{rmEffect(nCond)}).Properties.VariableNames=['ID' condNamesVerif{indEffect(1)} condNamesVerif{indEffect(2)} cond4effect{rmEffect(nCond)}];
    elseif numel(indEffect)==3
        tablemeans.(condNamesVerif{rmEffect(nCond)}).Properties.VariableNames=['ID' condNamesVerif{indEffect(1)} condNamesVerif{indEffect(2)} condNamesVerif{indEffect(3)} cond4effect{rmEffect(nCond)}];
    end
    clear tempTableMeans tempTable


    % Effect sizes
    phFieldnames=fieldnames(postHoc);
    % simple rm effects
    for o=1:size(order4ES.(condNamesVerif{rmEffect(nCond)}),1)
        dES{1}=means(:,order4ES.(condNamesVerif{rmEffect(nCond)})(o,1));
        dES{2}=means(:,order4ES.(condNamesVerif{rmEffect(nCond)})(o,2));
        [ES(o,1), ESsd(o,1), diff(o,1), diffSD(o,1), relDiff(o,1), relDiffSD(o,1)]=esDiffCalculation0D(dES,0);
    end
    postHoc.(condNamesVerif{rmEffect(nCond)})=[postHoc.(condNamesVerif{rmEffect(nCond)}) table(ES) table(diff) table(relDiff) table(ESsd) table(diffSD) table(relDiffSD)];
    clear ES diff relDiff ESsd diffSD relDiffSD

    % interaction of rm with independant effect
    if numel(indEffect)>0
        for nInd=1:numel(indEffect)
            [~,b]=sort(allModalities{indEffect(nInd)});
            for nMod=1:numel(modalitiesInd{nInd})
                for o=1:size(order4ES.(condNamesVerif{rmEffect(nCond)}),1)
                    dES{1}=means(idxIndependantEffect(:,nInd)==nMod,order4ES.(condNamesVerif{rmEffect(nCond)})(o,1));
                    dES{2}=means(idxIndependantEffect(:,nInd)==nMod,order4ES.(condNamesVerif{rmEffect(nCond)})(o,2));
                    [ES(o,b(nMod)), ESsd(o,b(nMod)), diff(o,b(nMod)), diffSD(o,b(nMod)), relDiff(o,b(nMod)), relDiffSD(o,b(nMod))]=esDiffCalculation0D(dES,0);
                end
            end
            ES=ES(:); diff=diff(:); relDiff=relDiff(:); ESsd=ESsd(:); diffSD=diffSD(:); relDiffSD=relDiffSD(:);
            postHoc.([condNamesVerif{rmEffect(nCond)} 'By' condNamesVerif{indEffect(nInd)}])=[postHoc.([condNamesVerif{rmEffect(nCond)} 'By' condNamesVerif{indEffect(nInd)}]) table(ES) table(diff) table(relDiff) table(ESsd) table(diffSD) table(relDiffSD)];
            clear ES diff relDiff ESsd diffSD relDiffSD
        end
    end

    % interaction of rm with douple independant
    if numel(indEffect)>1
        fieldNames=fieldnames(modInt2);
        for nInd=1:numel(fieldNames)
            [~,b]=sort(modInt2.(fieldNames{nInd}));
            for nMod=1:max(b)
                for o=1:size(order4ES.(condNamesVerif{rmEffect(nCond)}),1)
                    dES{1}=means(indInt2.(fieldNames{nInd})==nMod,order4ES.(condNamesVerif{rmEffect(nCond)})(o,1));
                    dES{2}=means(indInt2.(fieldNames{nInd})==nMod,order4ES.(condNamesVerif{rmEffect(nCond)})(o,2));
                    [ES(o,b(nMod)), ESsd(o,b(nMod)), diff(o,b(nMod)), diffSD(o,b(nMod)), relDiff(o,b(nMod)), relDiffSD(o,b(nMod))]=esDiffCalculation0D(dES,0);
                end
            end
            ES=ES(:); diff=diff(:); relDiff=relDiff(:); ESsd=ESsd(:); diffSD=diffSD(:); relDiffSD=relDiffSD(:);
            if ~isempty(findcolExact(phFieldnames,[condNamesVerif{rmEffect(nCond)} 'By' fieldNames{nInd}]))
                postHoc.([condNamesVerif{rmEffect(nCond)} 'By' fieldNames{nInd}])=[postHoc.([condNamesVerif{rmEffect(nCond)} 'By' fieldNames{nInd}]) table(ES) table(diff) table(relDiff) table(ESsd) table(diffSD) table(relDiffSD)];
            end
            clear ES diff relDiff ESsd diffSD relDiffSD
        end
    end
end

% interaction betwween rm effects
if numel(rmEffect)>1
    for nCond1=1:sum(isRM)
        for nCond2=1:sum(isRM)
            if nCond1~=nCond2
                for o=1:size(order4ES.([condNamesVerif{rmEffect(nCond1)} 'By'  condNamesVerif{rmEffect(nCond2)}]),1)
                    dES{1}=tXl(:,order4ES.([condNamesVerif{rmEffect(nCond1)} 'By'  condNamesVerif{rmEffect(nCond2)}])(o,1));
                    dES{2}=tXl(:,order4ES.([condNamesVerif{rmEffect(nCond1)} 'By'  condNamesVerif{rmEffect(nCond2)}])(o,2));
                    [ES(o,1), ESsd(o,1), diff(o,1), diffSD(o,1), relDiff(o,1), relDiffSD(o,1)]=esDiffCalculation0D(dES,0);
                end
                postHoc.([condNamesVerif{rmEffect(nCond1)} 'By'  condNamesVerif{rmEffect(nCond2)}])=[postHoc.([condNamesVerif{rmEffect(nCond1)} 'By'  condNamesVerif{rmEffect(nCond2)}]) table(ES) table(diff) table(relDiff) table(ESsd) table(diffSD) table(relDiffSD)];
                clear ES diff relDiff ESsd diffSD relDiffSD
            end
        end
    end
end

% interaction 2rm effects and 1IND
if numel(rmEffect)>1 & numel(indEffect)>0
    for nCond1=1:sum(isRM)
        for nCond2=1:sum(isRM)
            if nCond1~=nCond2
                for nInd=1:numel(indEffect)
                    for nModInd=1:numel(allModalities{indEffect(nInd)})
                        for o=1:size(order4ES.([condNamesVerif{rmEffect(nCond1)} 'By'  condNamesVerif{rmEffect(nCond2)}]),1)
                            dES{1}=tXl(idxIndependantEffect(:,nInd)==nModInd,order4ES.([condNamesVerif{rmEffect(nCond1)} 'By'  condNamesVerif{rmEffect(nCond2)}])(o,1));
                            dES{2}=tXl(idxIndependantEffect(:,nInd)==nModInd,order4ES.([condNamesVerif{rmEffect(nCond1)} 'By'  condNamesVerif{rmEffect(nCond2)}])(o,2));
                            [ES(o,nModInd), ESsd(o,nModInd), diff(o,nModInd), diffSD(o,nModInd), relDiff(o,nModInd), relDiffSD(o,nModInd)]=esDiffCalculation0D(dES,0);
                        end
                    end
                    ES=ES(:); diff=diff(:); relDiff=relDiff(:); ESsd=ESsd(:); diffSD=diffSD(:); relDiffSD=relDiffSD(:);
                    postHoc.([condNamesVerif{rmEffect(nCond1)} 'By'  condNamesVerif{rmEffect(nCond2)} 'By' condNamesVerif{indEffect(nInd)}])=[postHoc.([condNamesVerif{rmEffect(nCond1)} 'By'  condNamesVerif{rmEffect(nCond2)} 'By' condNamesVerif{indEffect(nInd)}]) table(ES) table(diff) table(relDiff) table(ESsd) table(diffSD) table(relDiffSD)];
                    clear ES diff relDiff ESsd diffSD relDiffSD
                end
            end
        end
    end
end


% simple independant effect
dataAllRm=nanmean(tXl,2);
if numel(indEffect)>0
    for nInd=1:numel(indEffect)
        for o=1:size(order4ES.(condNamesVerif{indEffect(nInd)}),1)
            dES{1}=dataAllRm(order4ES.(condNamesVerif{indEffect(nInd)})(o,1)==idxIndependantEffect(:,nInd));
            dES{2}=dataAllRm(order4ES.(condNamesVerif{indEffect(nInd)})(o,2)==idxIndependantEffect(:,nInd));
            [ES(o,1), ESsd(o,1), diff(o,1), diffSD(o,1), relDiff(o,1), relDiffSD(o,1)]=esDiffCalculation0D(dES,1);
        end
        postHoc.(condNamesVerif{indEffect(nInd)})=[postHoc.(condNamesVerif{indEffect(nInd)}) table(ES) table(diff) table(relDiff) table(ESsd) table(diffSD) table(relDiffSD)];
        clear ES diff relDiff ESsd diffSD relDiffSD
    end
end

% interaction of indpendant with rm effect
if numel(indEffect)>0 & numel(rmEffect)>0
    for nInd=1:numel(indEffect)
        for nRm=1:numel(rmEffect)
            for o=1:size(order4ES.(condNamesVerif{indEffect(nInd)}),1)
                for oRm=1:size(col4means{nRm},1)
                    dES{1}=nanmean(tXl(order4ES.(condNamesVerif{indEffect(nInd)})(o,1)==idxIndependantEffect(:,nInd),col4means{nRm}(oRm,:)),2);
                    dES{2}=nanmean(tXl(order4ES.(condNamesVerif{indEffect(nInd)})(o,2)==idxIndependantEffect(:,nInd),col4means{nRm}(oRm,:)),2);
                    [ES(o,oRm), ESsd(o,oRm), diff(o,oRm), diffSD(o,oRm), relDiff(o,oRm), relDiffSD(o,oRm)]=esDiffCalculation0D(dES,1);
                end
            end
            ES=ES(:); diff=diff(:); relDiff=relDiff(:); ESsd=ESsd(:); diffSD=diffSD(:); relDiffSD=relDiffSD(:);
            postHoc.([condNamesVerif{indEffect(nInd)} 'By' condNamesVerif{rmEffect(nRm)}])=[postHoc.([condNamesVerif{indEffect(nInd)} 'By' condNamesVerif{rmEffect(nRm)}]) table(ES) table(diff) table(relDiff) table(ESsd) table(diffSD) table(relDiffSD)];
            clear ES diff relDiff ESsd diffSD relDiffSD
        end
    end
end

% interaction of indpendant with douple rm
if numel(rmEffect)==2 & numel(indEffect)>0
    for nInd=1:numel(indEffect)
        for nRm1=1:numel(rmEffect)
            for nRm2=1:numel(rmEffect)
                if nRm2~=nRm1
                    %                     cols=unique([col4means{nRm1}(:) col4means{nRm2}(:)]);
                    % thinking about 3rd rm effect by averaging it
                    for o=1:size(order4ES.(condNamesVerif{indEffect(nInd)}),1)
                        for oRm=1:size(tXl,2)
                            dES{1}=(tXl(order4ES.(condNamesVerif{indEffect(nInd)})(o,1)==idxIndependantEffect(:,nInd),oRm));
                            dES{2}=(tXl(order4ES.(condNamesVerif{indEffect(nInd)})(o,2)==idxIndependantEffect(:,nInd),oRm));
                            [ES(o,oRm), ESsd(o,oRm), diff(o,oRm), diffSD(o,oRm), relDiff(o,oRm), relDiffSD(o,oRm)]=esDiffCalculation0D(dES,1);
                        end
                    end
                    ES=ES(:); diff=diff(:); relDiff=relDiff(:); ESsd=ESsd(:); diffSD=diffSD(:); relDiffSD=relDiffSD(:);
                    postHoc.([condNamesVerif{indEffect(nInd)} 'By' condNamesVerif{rmEffect(nRm1)} 'By' condNamesVerif{rmEffect(nRm2)}])=[postHoc.([condNamesVerif{indEffect(nInd)} 'By' condNamesVerif{rmEffect(nRm1)} 'By' condNamesVerif{rmEffect(nRm2)}]) table(ES) table(diff) table(relDiff) table(ESsd) table(diffSD) table(relDiffSD)];
                    clear ES diff relDiff ESsd diffSD relDiffSD
                end
            end
        end
    end
end

% interaction of double indpendant
if numel(indEffect)>1
    data4ES=nanmean(tXl,2);
    for nInd1=1:numel(indEffect)
        for nInd2=1:numel(indEffect)
            if nInd2~=nInd1
                [~,b]=sort(modalitiesInd{nInd2});
                for nMod2=1:numel(modalitiesInd{nInd2})
                    for o=1:size(order4ES.(condNamesVerif{indEffect(nInd1)}),1)
                        dES{1}=data4ES(order4ES.(condNamesVerif{indEffect(nInd1)})(o,1)==idxIndependantEffect(:,nInd1) & idxIndependantEffect(:,nInd2)==nMod2);
                        dES{2}=data4ES(order4ES.(condNamesVerif{indEffect(nInd1)})(o,2)==idxIndependantEffect(:,nInd1) & idxIndependantEffect(:,nInd2)==nMod2);
                        [ES(o,b(nMod2)), ESsd(o,b(nMod2)), diff(o,b(nMod2)), diffSD(o,b(nMod2)), relDiff(o,b(nMod2)), relDiffSD(o,b(nMod2))]=esDiffCalculation0D(dES,1);
                    end
                end
                ES=ES(:); diff=diff(:); relDiff=relDiff(:); ESsd=ESsd(:); diffSD=diffSD(:); relDiffSD=relDiffSD(:);
                postHoc.([condNamesVerif{indEffect(nInd1)} 'By' condNamesVerif{indEffect(nInd2)}])=[postHoc.([condNamesVerif{indEffect(nInd1)} 'By' condNamesVerif{indEffect(nInd2)}]) table(ES) table(diff) table(relDiff) table(ESsd) table(diffSD) table(relDiffSD)];
                clear ES diff relDiff ESsd diffSD relDiffSD
            end
        end
    end
end

% interaction of double indpendant with simple rm
if numel(indEffect)>1 & numel(rmEffect)>0
    for nRm=1:numel(rmEffect)
        for nInd1=1:numel(indEffect)
            for nInd2=1:numel(indEffect)
                if nInd1~=nInd2
                    for nModRm=1:size(col4means{nRm},1)
                        data4ES=nanmean(tXl(:,col4means{nRm}(nModRm,:)),2);
                        [~,bRm]=sort(modalitiesRM{nRm});
                        for nMod2=1:numel(modalitiesInd{nInd2})
                            for o=1:size(order4ES.(condNamesVerif{nInd1}),1)
                                dES{1}=data4ES(order4ES.(condNamesVerif{indEffect(nInd1)})(o,1)==idxIndependantEffect(:,nInd1) & idxIndependantEffect(:,nInd2)==nMod2);
                                dES{2}=data4ES(order4ES.(condNamesVerif{indEffect(nInd1)})(o,2)==idxIndependantEffect(:,nInd1) & idxIndependantEffect(:,nInd2)==nMod2);
                                [ES(o,nMod2,bRm(nModRm)), ESsd(o,nMod2,bRm(nModRm)), diff(o,nMod2,bRm(nModRm)), diffSD(o,nMod2,bRm(nModRm)), relDiff(o,nMod2,bRm(nModRm)), relDiffSD(o,nMod2,bRm(nModRm))]=esDiffCalculation0D(dES,1);
                            end
                        end
                    end
                    ES=ES(:); diff=diff(:); relDiff=relDiff(:); ESsd=ESsd(:); diffSD=diffSD(:); relDiffSD=relDiffSD(:);
                    if ~isempty(findcolExact(phFieldnames,[condNamesVerif{nInd1} 'By' condNamesVerif{nInd2} 'By' condNamesVerif{rmEffect(nRm)}]))
                        postHoc.([condNamesVerif{nInd1} 'By' condNamesVerif{nInd2} 'By' condNamesVerif{rmEffect(nRm)}])=[postHoc.([condNamesVerif{nInd1} 'By' condNamesVerif{nInd2} 'By' condNamesVerif{rmEffect(nRm)}]) table(ES) table(diff) table(relDiff) table(ESsd) table(diffSD) table(relDiffSD)];
                    end
                    clear ES diff relDiff ESsd diffSD relDiffSD
                end
            end
        end
    end
end

% to do 2 for ES
% interaction of double indpendant with double rm
% all interactions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SAVING TABLES
fieldNames4SavePH=fieldnames(postHoc);
fieldNames4SaveMeans=fieldnames(tablemeans);

tablemeans=AOVreformTable(tablemeans, stats, allMod, aov);

Tables.means=tablemeans;
Tables.aov=aov;
Tables.posthoc=postHoc;


if ~isempty(saveDir)
    writetable(aov,fullfile(saveDir, 'Tables.xlsx'),"Sheet",'AOV','WriteMode', 'replacefile')
    for i=1:numel(fieldNames4SaveMeans)
        writetable(tablemeans.(fieldNames4SaveMeans{i}),fullfile(saveDir, 'Tables.xlsx'),"Sheet",['MEANS_' fieldNames4SaveMeans{i}])
    end
    for i=1:numel(fieldNames4SavePH)
        writetable(postHoc.(fieldNames4SavePH{i}),fullfile(saveDir, 'Tables.xlsx'),"Sheet",['PH_' fieldNames4SavePH{i}])
    end
    save(fullfile(saveDir, 'Tables'),'Tables')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% PLOT
if ~isempty(saveDir)

    if ~isempty(findcol(fieldnames(stats),"addUnivariate"))
        for nRm=1:numel(modalitiesRM)
            modalitiesRM{nRm}(end)=[];
            allMod.(condNamesVerif{rmEffect(nRm)})(end)=[];
            maxN=numel(modalitiesRM{nRm})+1;
            for n=1:size(order4ES.(condNamesVerif{rmEffect(nRm)}),1)
                ESmax=max(order4ES.(condNamesVerif{rmEffect(nRm)})(n,:));
                if ESmax==maxN
                    order4ES.(condNamesVerif{rmEffect(nRm)})(n,2)=order4ES.(condNamesVerif{rmEffect(nRm)})(n,1);
                end
            end
        end
    end

    for nCond=1:numel(condNames)
        if statsLines
            mkdir(fullfile(saveDir, 'Lines', condNames{nCond}))
        end
        mkdir(fullfile(saveDir, 'Text', condNames{nCond}))
    end

    %% Repeated measures effects
    %% MAIN EFFECT
    % Lines
    if ~isempty(rmEffect)
        if statsLines
            for nRm=1:numel(modalitiesRM)

                f=figure('units','centimeters','position',[0 0 6+4*numel(modalitiesRM{nRm}) 4+9/16*4*numel(modalitiesRM{nRm})],'visible','off');

                for x=1:numel(modalitiesRM{nRm})
                    dataMeans(x)=nanmean(nanmean(data4plot.allData(:,col4means{nRm}(x,:)),2));
                    dataSD(x)=nanstd(nanmean(data4plot.allData(:,col4means{nRm}(x,:)),2));
                end

                for x=1:numel(modalitiesRM{nRm})

                    h=bar(x,dataMeans(x)); hold on
                    h(1).FaceColor=colors{rmEffect(nRm)}(x,:);

                    if plotSD==1
                        if dataMeans(x)>=0
                            SD(1)=0;
                            SD(2)=dataSD(x);
                        else
                            SD(1)=dataSD(x);
                            SD(2)=0;
                        end
                        errorbar(x,dataMeans(x),SD(1),SD(2),'k','LineStyle','none')
                    end
                end

                if indivLines==1
                    xl=1:numel(modalitiesRM{nRm});
                    for xx=1:numel(modalitiesRM{nRm})
                        dataLines(:,xx)=nanmean(data4plot.allData(:,col4means{nRm}(xx,:)),2);
                    end
                    plot(xl,dataLines,'k--'); hold on
                    scatter(xl,dataLines,'k+')
                end
                clear dataLines

                for x=1:numel(modalitiesRM{nRm})
                    if abs(dataMeans(x))<1
                        text(x,0.5*(dataMeans(x)),sprintf('%.3f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                    elseif abs(dataMeans(x))<10
                        text(x,0.5*(dataMeans(x)),sprintf('%.2f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                    elseif abs(dataMeans(x))<100
                        text(x,0.5*(dataMeans(x)),sprintf('%.1f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                    else
                        text(x,0.5*(dataMeans(x)),sprintf('%.0f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                    end
                end

                xticks(1:numel(modalitiesRM{nRm}))
                xticklabels(modalitiesRM{nRm})
                xlabel(condNames{rmEffect(nRm)})
                ylabel(units)
                yl=ylim;
                box off
                ax=gca;
                ax.XGrid='off';
                ax.YGrid='on';

                nAov=findcolExact(aov.Effect, verifFieldName(condNames{rmEffect(nRm)}));
                pMAIN=aov{nAov,6};

                amp=[max(yl)-min(yl)];
                isSignificant=0;
                pValues=ones(1,size(postHoc.(condNamesVerif{rmEffect(nRm)}),1));
                pSelected=0;
                if pMAIN<pcritical(1) &  pSelected==0
                    isSignificant=1;
                    pValues=postHoc.(condNamesVerif{rmEffect(nRm)}){:,5};
                    pSelected=1;
                end

                if pMAIN<pcritical(end) &  pSelected==0
                    pValues=postHoc.(condNamesVerif{rmEffect(nRm)}){:,5};
                    pSelected=1;
                end

                if pSelected==1
                    nSignificant=1;
                    for i=1:numel(pValues)
                        pV=pValues(i);
                        if pV<=pcritical(end)
                            if abs(yl(2))>abs(yl(1))
                                if pV<pcritical(1) & isSignificant==1
                                    hline(yl(2)+0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNamesVerif{rmEffect(nRm)})(i,:),'lineWidth',1.5);
                                else
                                    hline(yl(2)+0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNamesVerif{rmEffect(nRm)})(i,:),'lineWidth',1.5);
                                end
                            else
                                if pV<pcritical(1) & isSignificant==1
                                    hline(yl(1)-0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNamesVerif{rmEffect(nRm)})(i,:),'lineWidth',1.5);
                                else
                                    hline(yl(1)-0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNamesVerif{rmEffect(nRm)})(i,:),'lineWidth',1.5);
                                end
                            end
                            nSignificant=nSignificant+1;
                        end
                    end
                end

                print('-dtiff',['-r' num2str(imageResolution)],fullfile(saveDir, 'Lines', condNames{rmEffect(nRm)}, 'All participants'))
                close

                if indivLines==2
                    for nInd=1:numel(indEffect)
                        f=figure('units','centimeters','position',[0 0 6+4*numel(modalitiesRM{nRm}) 4+9/16*1.5*numel(modalitiesRM{nRm})],'visible','off');
                        for nMod=1:numel(allMod.(condNamesVerif{indEffect(nInd)}))
                            for xl=1:numel(modalitiesRM{nRm})
                                selectIndivLines(:,xl)=mean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(:,col4means{nRm}(xl,:)),2);
                            end
                            xl=1:numel(modalitiesRM{nRm});
                            for n=1:size(selectIndivLines,1)
                                if n==1
                                    plot(xl,selectIndivLines(n,:),'color', colors{indEffect(nInd)}(nMod,:),'linewidth',1.05,'linestyle','-'); hold on
                                    scatter(xl,selectIndivLines(n,:), '+', 'MarkerFaceColor', colors{indEffect(nInd)}(nMod,:), 'MarkerEdgeColor', colors{indEffect(nInd)}(nMod,:),'handlevisibility','off')
                                else
                                    plot(xl,selectIndivLines(n,:),'color', colors{indEffect(nInd)}(nMod,:),'linewidth',1.05,'linestyle','-', 'handlevisibility','off'); hold on
                                    scatter(xl,selectIndivLines(n,:), '+', 'MarkerFaceColor', colors{indEffect(nInd)}(nMod,:), 'MarkerEdgeColor', colors{indEffect(nInd)}(nMod,:),'handlevisibility','off')
                                end
                            end
                            clear selectIndivLines
                        end
                        xlim(xlp)
                        xticks(1:numel(modalitiesRM{nRm}))
                        xticklabels(modalitiesRM{nRm})
                        xlabel(condNames{rmEffect(nRm)})
                        ylabel(units)
                        yl=ylim;
                        box off
                        ax=gca;
                        ax.XGrid='off';
                        ax.YGrid='on';
                        legend(allModalities{indEffect(nInd)}, "box", "off", "Location", "northeast");

                        print('-dtiff',['-r' num2str(imageResolution)],fullfile(saveDir, 'Lines', condNames{rmEffect(nRm)}, ['All Participants and ' condNames{indEffect(nInd)}]))
                        close
                    end
                end

                clear ax yl dataMeans dataSD
            end
        end

        % Text
        for nRm=1:numel(effectRM)
            f=figure('units','centimeters','position',[0 0 6+4*numel(modalitiesRM{nRm}) 4+9/16*4*numel(modalitiesRM{nRm})],'visible','off');

            for x=1:numel(modalitiesRM{nRm})
                dataMeans(x)=nanmean(nanmean(data4plot.allData(:,col4means{nRm}(x,:)),2));
                dataSD(x)=nanstd(nanmean(data4plot.allData(:,col4means{nRm}(x,:)),2));
            end

            for x=1:numel(modalitiesRM{nRm})

                h=bar(x,dataMeans(x)); hold on
                h(1).FaceColor=colors{rmEffect(nRm)}(x,:);

                if plotSD==1
                    if dataMeans(x)>=0
                        SD(1)=0;
                        SD(2)=dataSD(x);
                    else
                        SD(1)=dataSD(x);
                        SD(2)=0;
                    end
                    errorbar(x,dataMeans(x),SD(1),SD(2),'k','LineStyle','none')
                end
            end

            if indivLines==1
                xl=1:numel(modalitiesRM{nRm});
                for xx=1:numel(modalitiesRM{nRm})
                    dataLines(:,xx)=nanmean(data4plot.allData(:,col4means{nRm}(xx,:)),2);
                end
                plot(xl,dataLines,'k--'); hold on
                scatter(xl,dataLines,'k+')
            end

            for x=1:numel(modalitiesRM{nRm})
                if abs(dataMeans(x))<1
                    text(x,0.5*(dataMeans(x)),sprintf('%.3f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                elseif abs(dataMeans(x))<10
                    text(x,0.5*(dataMeans(x)),sprintf('%.2f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                elseif abs(dataMeans(x))<100
                    text(x,0.5*(dataMeans(x)),sprintf('%.1f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                else
                    text(x,0.5*(dataMeans(x)),sprintf('%.0f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                end
            end

            xticks(1:numel(modalitiesRM{nRm}))
            xticklabels(modalitiesRM{nRm})
            xlabel(condNames{rmEffect(nRm)})
            ylabel(units)
            yl=ylim;
            box off
            ax=gca;
            ax.XGrid='off';
            ax.YGrid='on';

            nAov=findcolExact(aov.Effect, verifFieldName(condNames{rmEffect(nRm)}));
            pMAIN=aov{nAov,6};

            amp=[max(yl)-min(yl)];
            isSignificant=0;
            pValues=ones(1,size(postHoc.(condNamesVerif{rmEffect(nRm)}),1));
            pSelected=0;
            if pMAIN<pcritical(1) &  pSelected==0
                isSignificant=1;
                pValues=postHoc.(condNamesVerif{rmEffect(nRm)}){:,5};
                pSelected=1;
            end
            if pMAIN<pcritical(end) &  pSelected==0
                pValues=postHoc.(condNamesVerif{rmEffect(nRm)}){:,5};
                pSelected=1;
            end

            dataMeans4pv=dataMeans;
            if plotSD==1
                dataMeans4pv=dataMeans+sign(dataMeans).*dataSD;
            end
            if indivLines==1
                dataMeans4pv=max(abs(dataLines)).*sign(dataMeans);
            end
            if plotSD==1 && indivLines==1
                dataMeans4pv=sign(dataMeans).*max([abs(dataMeans+sign(dataMeans).*dataSD); max(abs(dataLines))]);
            end
            clear dataLines

            nColSignificant=ones(1,numel(dataMeans));
            nSignificant=1;
            if pSelected==1
                for i=1:numel(pValues)
                    pV=pValues(i);
                    if pV<=pcritical(end)
                        if abs(yl(2))>abs(yl(1))
                            if abs(dataMeans4pv(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1)))>abs(dataMeans4pv(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2)))
                                addPvalue(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1), dataMeans4pv(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1))+0.055*nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1))*amp, pV, pcritical, colors{rmEffect(nRm)}(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2),:))
                                nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1))=nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1))+1;
                            else
                                addPvalue(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2), dataMeans4pv(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2))+0.055*nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2))*amp, pV, pcritical, colors{rmEffect(nRm)}(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1),:))
                                nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2))=nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2))+1;
                            end
                        else
                            if abs(dataMeans4pv(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1)))>abs(dataMeans4pv(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2)))
                                addPvalue(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1), dataMeans4pv(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1))-0.055*nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1))*amp, pV, pcritical, colors{rmEffect(nRm)}(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2),:))
                                nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1))=nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1))+1;
                            else
                                addPvalue(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2), dataMeans4pv(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2))-0.055*nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2))*amp, pV, pcritical, colors{rmEffect(nRm)}(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1),:))
                                nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2))=nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2))+1;
                            end
                        end
                        nSignificant=nSignificant+1;
                    end
                end
            end

            yl=yl*1.05.^max(nColSignificant);
            set(ax,'ylim',[min(yl) max(yl)])

            xlp=xlim;

            print('-dtiff',['-r' num2str(imageResolution)],fullfile(saveDir, 'Text', condNames{rmEffect(nRm)}, 'All participants'))
            close

            if indivLines==2
                for nInd=1:numel(indEffect)
                    f=figure('units','centimeters','position',[0 0 6+4*numel(modalitiesRM{nRm}) 4+9/16*1.5*numel(modalitiesRM{nRm})],'visible','off');
                    for nMod=1:numel(allMod.(condNamesVerif{indEffect(nInd)}))
                        for xl=1:numel(modalitiesRM{nRm})
                            selectIndivLines(:,xl)=mean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(:,col4means{nRm}(xl,:)),2);
                        end
                        hline(0, 'xLimits',[0 (numel(modalitiesRM{nRm})+1)],'linetype','-k', 'linewidth',0.5)
                        xl=1:numel(modalitiesRM{nRm});
                        for n=1:size(selectIndivLines,1)
                            if n==1
                                plot(xl,selectIndivLines(n,:),'color', colors{indEffect(nInd)}(nMod,:),'linewidth',1.05,'linestyle','-'); hold on
                                scatter(xl,selectIndivLines(n,:), '+', 'MarkerFaceColor', colors{indEffect(nInd)}(nMod,:), 'MarkerEdgeColor', colors{indEffect(nInd)}(nMod,:),'handlevisibility','off')
                            else
                                plot(xl,selectIndivLines(n,:),'color', colors{indEffect(nInd)}(nMod,:),'linewidth',1.05,'linestyle','-', 'handlevisibility','off'); hold on
                                scatter(xl,selectIndivLines(n,:), '+', 'MarkerFaceColor', colors{indEffect(nInd)}(nMod,:), 'MarkerEdgeColor', colors{indEffect(nInd)}(nMod,:),'handlevisibility','off')
                            end
                        end
                        clear selectIndivLines
                    end
                    xlim(xlp)
                    xticks(1:numel(modalitiesRM{nRm}))
                    xticklabels(modalitiesRM{nRm})
                    xlabel(condNames{rmEffect(nRm)})
                    ylabel(units)
                    yl=ylim;
                    ylim([-max(abs(yl)) max(abs(yl))])
                    box off
                    ax=gca;
                    ax.XGrid='off';
                    ax.YGrid='on';
                    legend(allModalities{indEffect(nInd)}, "box", "off", "Location", "northeast");

                    print('-dtiff',['-r' num2str(imageResolution)],fullfile(saveDir, 'Text', condNames{rmEffect(nRm)}, ['All Participants and ' condNames{indEffect(nInd)}]))
                    close
                end
            end

            clear ax yl dataMeans dataSD

        end
    end

    %% INTERACTIONS of RM
    % lines
    if ~isempty(rmEffect)
        if numel(effectRM)>1
            if statsLines
                for nRm1=1:numel(modalitiesRM)
                    for nRm2=1:numel(modalitiesRM)
                        if nRm1~=nRm2

                            f=figure('units','centimeters','position',[0 0 6+4*numel(allModalities{rmEffect(nRm1)}) 4+numel(cond4effect{rmEffect(nRm2)})*numel(cond4effect{rmEffect(nRm1)})*9/16*4],'visible','off');

                            for nModRm=1:numel(cond4effect{rmEffect(nRm2)})

                                varData=tXl(:,col4means{nRm1}(:,nModRm));
                                dataMeans=nanmean(varData);
                                dataSD=nanstd(varData);

                                subplot(numel(cond4effect{rmEffect(nRm2)}),1,nModRm);
                                for x=1:numel(cond4effect{rmEffect(nRm1)})
                                    h=bar(x,dataMeans(x)); hold on
                                    h(1).FaceColor=colors{rmEffect(nRm1)}(x,:);

                                    if plotSD==1
                                        if dataMeans(x)>=0
                                            SD(1)=0;
                                            SD(2)=dataSD(x);
                                        else
                                            SD(1)=dataSD(x);
                                            SD(2)=0;
                                        end
                                        errorbar(x,dataMeans(x),SD(1),SD(2),'k','LineStyle','none')
                                    end
                                end
                                if indivLines==1
                                    xl=1:numel(modalitiesRM{nRm1});
                                    plot(xl,varData,'k--'); hold on
                                    scatter(xl,varData,'k+')
                                end
                                for x=1:numel(cond4effect{rmEffect(nRm1)})
                                    if abs(dataMeans(x))<1
                                        text(x,0.5*(dataMeans(x)),sprintf('%.3f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                    elseif abs(dataMeans(x))<10
                                        text(x,0.5*(dataMeans(x)),sprintf('%.2f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                    elseif abs(dataMeans(x))<100
                                        text(x,0.5*(dataMeans(x)),sprintf('%.1f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                    else
                                        text(x,0.5*(dataMeans(x)),sprintf('%.0f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                    end
                                end

                                xticks(1:numel(cond4effect{rmEffect(nRm1)}))
                                xticklabels(cond4effect{rmEffect(nRm1)})
                                if nModRm==numel(cond4effect{rmEffect(nRm2)})
                                    xlabel(condNames{rmEffect(nRm1)})
                                end
                                ylabel(units)
                                box off
                                ax{nModRm}=gca;
                                ax{nModRm}.XGrid='off';
                                ax{nModRm}.YGrid='on';
                                yl(:,nModRm)=ylim;
                                title(cond4effect{rmEffect(nRm2)}{nModRm})

                            end

                            nAovInt=findcolExact(aov.Effect,[condNamesVerif{rmEffect(nRm1)} ':' condNamesVerif{rmEffect(nRm2)}]);
                            if isempty(nAovInt)
                                nAovInt=findcolExact(aov.Effect,[condNamesVerif{rmEffect(nRm2)} ':' condNamesVerif{rmEffect(nRm1)}]);
                            end
                            nAov=findcolExact(aov.Effect, verifFieldName(condNamesVerif{rmEffect(nRm1)}));
                            pINT=aov{nAovInt,6};
                            pMAIN=aov{nAov,6};

                            amp=[max(max(yl))-min(min(yl))];
                            for nModRm=1:numel(cond4effect{rmEffect(nRm2)})

                                varData=tXl(:,col4means{nRm1}(:,nModRm));

                                isSignificant=0;
                                set(f,'CurrentAxes',ax{nModRm});
                                pValues=ones(1,size(postHoc.(condNamesVerif{rmEffect(nRm1)}),1));
                                pSelected=0;
                                if pINT<pcritical(1)
                                    isSignificant=1;
                                    phCut=findcolExact(postHoc.([condNamesVerif{rmEffect(nRm1)} 'By' condNamesVerif{rmEffect(nRm2)}]){:,1},allMod.(condNamesVerif{rmEffect(nRm2)}){nModRm});
                                    pValues=postHoc.([condNamesVerif{rmEffect(nRm1)} 'By' condNamesVerif{rmEffect(nRm2)}]){phCut,6};
                                    pSelected=1;
                                end
                                %                         if pMAIN<pcritical(1) &  pSelected==0
                                %                             isSignificant=1;
                                %                             pValues=postHoc.(condNamesVerif{rmEffect(nRm1)}){:,5};
                                %                             pSelected=1;
                                %                         end
                                if pINT<pcritical(end) &  pSelected==0
                                    phCut=findcolExact(postHoc.([condNamesVerif{rmEffect(nRm1)} 'By' condNamesVerif{rmEffect(nRm2)}]){:,1},allMod.(condNamesVerif{rmEffect(nRm2)}){nModRm});
                                    pValues=postHoc.([condNamesVerif{rmEffect(nRm1)} 'By' condNamesVerif{rmEffect(nRm2)}]){phCut,6};
                                    pSelected=1;
                                end
                                %                         if pMAIN<pcritical(end) &  pSelected==0
                                %                             pValues=postHoc.(condNamesVerif{rmEffect(nRm1)}){:,5};
                                %                             pSelected=1;
                                %                         end

                                if pSelected==1
                                    nSignificant=1;
                                    for i=1:numel(pValues)
                                        pV=pValues(i);
                                        if pV<=pcritical(end)
                                            if yl(2,nModRm)>0
                                                if pV<pcritical(1) & isSignificant==1
                                                    hline(yl(2,nModRm)+0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNamesVerif{rmEffect(nRm1)})(i,:),'lineWidth',1.5);
                                                else
                                                    hline(yl(2,nModRm)+0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNamesVerif{rmEffect(nRm1)})(i,:),'lineWidth',1.5);
                                                end
                                            else
                                                if pV<pcritical(1) & isSignificant==1
                                                    hline(yl(1,nModRm)-0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNamesVerif{rmEffect(nRm1)})(i,:),'lineWidth',1.5);
                                                else
                                                    hline(yl(1,nModRm)-0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNamesVerif{rmEffect(nRm1)})(i,:),'lineWidth',1.5);
                                                end
                                            end
                                            nSignificant=nSignificant+1;
                                        end
                                    end
                                end

                                yl(:,nModRm)=ylim;

                            end

                            for nModRm=1:numel(ax)
                                set(ax{nModRm},'ylim',1.05*[min(min(yl)) max(max(yl))])
                            end

                            print('-dtiff',['-r' num2str(imageResolution)],fullfile(saveDir, 'Lines', condNames{rmEffect(nRm1)}, ['All Participants by ' condNames{rmEffect(nRm2)}]))
                            close
                            clear ax yl dataMeans dataSD

                        end
                    end
                end
            end

            % text
            for nRm1=1:numel(modalitiesRM)
                for nRm2=1:numel(modalitiesRM)
                    if nRm1~=nRm2

                        f=figure('units','centimeters','position',[0 0 6+4*numel(allModalities{rmEffect(nRm1)}) 4+numel(cond4effect{rmEffect(nRm2)})*numel(cond4effect{rmEffect(nRm1)})*9/16*4],'visible','off');
                        for nModRm=1:numel(cond4effect{rmEffect(nRm2)})

                            varData=tXl(:,col4means{nRm1}(:,nModRm));
                            dataMeans=nanmean(varData);
                            dataSD=nanstd(varData);

                            subplot(numel(cond4effect{rmEffect(nRm2)}),1,nModRm);
                            for x=1:numel(cond4effect{rmEffect(nRm1)})
                                h=bar(x,dataMeans(x)); hold on
                                h(1).FaceColor=colors{rmEffect(nRm1)}(x,:);

                                if plotSD==1
                                    if dataMeans(x)>=0
                                        SD(1)=0;
                                        SD(2)=dataSD(x);
                                    else
                                        SD(1)=dataSD(x);
                                        SD(2)=0;
                                    end
                                    errorbar(x,dataMeans(x),SD(1),SD(2),'k','LineStyle','none')
                                end
                            end

                            if indivLines==1
                                xl=1:numel(modalitiesRM{nRm1});
                                plot(xl,varData,'k--'); hold on
                                scatter(xl,varData,'k+')
                            end

                            for x=1:numel(cond4effect{rmEffect(nRm1)})
                                if abs(dataMeans(x))<1
                                    text(x,0.5*(dataMeans(x)),sprintf('%.3f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                elseif abs(dataMeans(x))<10
                                    text(x,0.5*(dataMeans(x)),sprintf('%.2f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                elseif abs(dataMeans(x))<100
                                    text(x,0.5*(dataMeans(x)),sprintf('%.1f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                else
                                    text(x,0.5*(dataMeans(x)),sprintf('%.0f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                end
                            end

                            xticks(1:numel(cond4effect{rmEffect(nRm1)}))
                            xticklabels(cond4effect{rmEffect(nRm1)})
                            if nModRm==numel(cond4effect{rmEffect(nRm2)})
                                xlabel(condNames{rmEffect(nRm1)})
                            end
                            ylabel(units)
                            box off
                            ax{nModRm}=gca;
                            ax{nModRm}.XGrid='off';
                            ax{nModRm}.YGrid='on';
                            yl(:,nModRm)=ylim;
                            title(cond4effect{rmEffect(nRm2)}{nModRm})

                        end

                        nAovInt=findcolExact(aov.Effect,[condNamesVerif{rmEffect(nRm1)} ':' condNamesVerif{rmEffect(nRm2)}]);
                        if isempty(nAovInt)
                            nAovInt=findcolExact(aov.Effect,[condNamesVerif{rmEffect(nRm2)} ':' condNamesVerif{rmEffect(nRm1)}]);
                        end
                        nAov=findcolExact(aov.Effect, verifFieldName(condNamesVerif{rmEffect(nRm1)}));
                        pINT=aov{nAovInt,6};
                        pMAIN=aov{nAov,6};

                        amp=[max(max(yl))-min(min(yl))];
                        for nModRm=1:numel(cond4effect{rmEffect(nRm2)})

                            varData=tXl(:,col4means{nRm1}(:,nModRm));
                            dataMeans=nanmean(varData);
                            dataSD=nanstd(varData);
                            dataMeans4pv=dataMeans;
                            if plotSD==1
                                dataMeans4pv=dataMeans+sign(dataMeans).*dataSD;
                            end
                            if indivLines==1
                                dataMeans4pv=max(abs(varData)).*sign(dataMeans);
                            end
                            if plotSD==1 && indivLines==1
                                dataMeans4pv=sign(dataMeans).*max([abs(dataMeans+sign(dataMeans).*dataSD); max(abs(varData)).*sign(dataMeans)]);
                            end

                            isSignificant=0;
                            set(f,'CurrentAxes',ax{nModRm});
                            pValues=ones(1,size(postHoc.(condNamesVerif{rmEffect(nRm1)}),1));
                            pSelected=0;
                            if pINT<pcritical(1)
                                isSignificant=1;
                                phCut=findcolExact(postHoc.([condNamesVerif{rmEffect(nRm1)} 'By' condNamesVerif{rmEffect(nRm2)}]){:,1},allMod.(condNamesVerif{rmEffect(nRm2)}){nModRm});
                                pValues=postHoc.([condNamesVerif{rmEffect(nRm1)} 'By' condNamesVerif{rmEffect(nRm2)}]){phCut,6};
                                pSelected=1;
                            end
                            %                     if pMAIN<pcritical(1) &  pSelected==0
                            %                         isSignificant=1;
                            %                         pValues=postHoc.(condNamesVerif{rmEffect(nRm1)}){:,5};
                            %                         pSelected=1;
                            %                     end
                            if pINT<pcritical(end) &  pSelected==0
                                phCut=findcolExact(postHoc.([condNamesVerif{rmEffect(nRm1)} 'By' condNamesVerif{rmEffect(nRm2)}]){:,1},allMod.(condNamesVerif{rmEffect(nRm2)}){nModRm});
                                pValues=postHoc.([condNamesVerif{rmEffect(nRm1)} 'By' condNamesVerif{rmEffect(nRm2)}]){phCut,6};
                                pSelected=1;
                            end
                            %                     if pMAIN<pcritical(end) &  pSelected==0
                            %                         pValues=postHoc.(condNamesVerif{rmEffect(nRm1)}){:,5};
                            %                         pSelected=1;
                            %                     end

                            nSignificant=1;
                            nColSignificant=ones(1,numel(dataMeans));
                            if pSelected==1
                                for i=1:numel(pValues)
                                    pV=pValues(i);
                                    if pV<=pcritical(end)
                                        if abs(yl(2))>abs(yl(1))
                                            if abs(dataMeans(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1)))>abs(dataMeans(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2)))
                                                addPvalue(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1), dataMeans4pv(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1))+0.0225*numel(cond4effect{rmEffect(nRm2)})*nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1))*amp, pV, pcritical, colors{rmEffect(nRm1)}(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2),:))
                                                nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1))=nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1))+1;
                                            else
                                                addPvalue(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2), dataMeans4pv(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2))+0.0225*numel(cond4effect{rmEffect(nRm2)})*nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2))*amp, pV, pcritical, colors{rmEffect(nRm1)}(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1),:))
                                                nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2))=nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2))+1;
                                            end
                                        else
                                            if abs(dataMeans(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1)))>abs(dataMeans(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2)))
                                                addPvalue(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1), dataMeans4pv(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1))-0.0225*numel(cond4effect{rmEffect(nRm2)})*nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1))*amp, pV, pcritical, colors{rmEffect(nRm1)}(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2),:))
                                                nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1))=nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1))+1;
                                            else
                                                addPvalue(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2), dataMeans4pv(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2))-0.0225*numel(cond4effect{rmEffect(nRm2)})*nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2))*amp, pV, pcritical, colors{rmEffect(nRm1)}(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1),:))
                                                nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2))=nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2))+1;
                                            end
                                        end
                                        nSignificant=nSignificant+1;
                                    end
                                end
                            end

                            yl(:,nModRm)=ylim*1.05.^max(nColSignificant);

                        end

                        for nModRm=1:numel(ax)
                            set(ax{nModRm},'ylim',[min(min(yl)) max(max(yl))])
                        end

                        print('-dtiff',['-r' num2str(imageResolution)],fullfile(saveDir, 'Text', condNames{rmEffect(nRm1)}, ['All Participants by ' condNames{rmEffect(nRm2)}]))
                        close
                        clear ax yl dataMeans dataSD

                    end
                end
            end
        end
    end

    %% Main RM with simple IND interaction
    % Lines
    if ~isempty(rmEffect) & numel(indEffect)>0
        if statsLines
            for nRm=1:numel(effectRM)
                for nInd=1:numel(indEffect)
                    for nMod=1:numel(allMod.(condNamesVerif{indEffect(nInd)}))

                        f=figure('units','centimeters','position',[0 0 6+4*numel(modalitiesRM{nRm}) 4+9/16*4*numel(modalitiesRM{nRm})],'visible','off');

                        for x=1:numel(modalitiesRM{nRm})
                            dataMeans(x)=nanmean(nanmean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(:,col4means{nRm}(x,:)),2));
                            dataSD(x)=nanstd(nanmean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(:,col4means{nRm}(x,:)),2));
                        end

                        for x=1:numel(modalitiesRM{nRm})

                            h=bar(x,dataMeans(x)); hold on
                            h(1).FaceColor=colors{rmEffect(nRm)}(x,:);

                            if plotSD==1
                                if dataMeans(x)>=0
                                    SD(1)=0;
                                    SD(2)=dataSD(x);
                                else
                                    SD(1)=dataSD(x);
                                    SD(2)=0;
                                end
                                errorbar(x,dataMeans(x),SD(1),SD(2),'k','LineStyle','none')
                            end
                        end

                        if indivLines==1
                            for xl=1:numel(modalitiesRM{nRm})
                                selectIndivLines(:,xl)=mean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(:,col4means{nRm}(xl,:)),2);
                            end
                            xl=1:numel(modalitiesRM{nRm});
                            plot(xl,selectIndivLines,'k--'); hold on
                            scatter(xl,selectIndivLines,'k+')
                        end
                        clear selectIndivLines

                        for x=1:numel(modalitiesRM{nRm})
                            if abs(dataMeans(x))<1
                                text(x,0.5*(dataMeans(x)),sprintf('%.3f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                            elseif abs(dataMeans(x))<10
                                text(x,0.5*(dataMeans(x)),sprintf('%.2f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                            elseif abs(dataMeans(x))<100
                                text(x,0.5*(dataMeans(x)),sprintf('%.1f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                            else
                                text(x,0.5*(dataMeans(x)),sprintf('%.0f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                            end

                        end
                        xticks(1:numel(modalitiesRM{nRm}))
                        xticklabels(modalitiesRM{nRm})
                        xlabel(condNames{rmEffect(nRm)})
                        ylabel(units)
                        yl=ylim;
                        box off
                        ax=gca;
                        ax.XGrid='off';
                        ax.YGrid='on';

                        nAovInt=findcolExact(aov.Effect,[condNamesVerif{indEffect(nInd)} ':' condNamesVerif{rmEffect(nRm)}]);
                        nAov=findcolExact(aov.Effect, verifFieldName(condNamesVerif{rmEffect(nRm)}));
                        pINT=aov{nAovInt,6};
                        pMAIN=aov{nAov,6};

                        amp=[max(yl)-min(yl)];
                        isSignificant=0;

                        pValues=ones(1,size(postHoc.(condNamesVerif{rmEffect(nRm)}),1));
                        pSelected=0;
                        if pINT<pcritical(1)
                            isSignificant=1;
                            rows=findcolExact(postHoc.([condNamesVerif{rmEffect(nRm)} 'By' condNamesVerif{indEffect(nInd)}]){:,1},allMod.(condNamesVerif{indEffect(nInd)}){nMod});
                            pValues=postHoc.([condNamesVerif{rmEffect(nRm)} 'By' condNamesVerif{indEffect(nInd)}]){rows,6};
                            pSelected=1;
                        end
                        %                     if pMAIN<pcritical(1) &  pSelected==0
                        %                         isSignificant=1;
                        %                         pValues=postHoc.(condNamesVerif{rmEffect(nRm)}){:,5};
                        %                         pSelected=1;
                        %                     end
                        if pINT<pcritical(end) &  pSelected==0
                            rows=findcolExact(postHoc.([condNamesVerif{rmEffect(nRm)} 'By' condNamesVerif{indEffect(nInd)}]){:,1},allMod.(condNamesVerif{indEffect(nInd)}){nMod});
                            pValues=postHoc.([condNamesVerif{rmEffect(nRm)} 'By' condNamesVerif{indEffect(nInd)}]){rows,6};
                            pSelected=1;
                        end
                        %                     if pMAIN<pcritical(end) &  pSelected==0
                        %                         pValues=postHoc.(condNamesVerif{rmEffect(nRm)}){:,5};
                        %                         pSelected=1;
                        %                     end

                        if plotSD==1
                            dataMeans=dataMeans+sign(dataMeans).*dataSD;
                        end

                        nSignificant=1;
                        if pSelected==1
                            nSignificant=1;
                            for i=1:numel(pValues)
                                pV=pValues(i);
                                if pV<=pcritical(end)
                                    if abs(yl(2))>abs(yl(1))
                                        if pV<pcritical(1) & isSignificant==1
                                            hline(yl(2)+0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNamesVerif{rmEffect(nRm)})(i,:),'lineWidth',1.5);
                                        else
                                            hline(yl(2)+0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNamesVerif{rmEffect(nRm)})(i,:),'lineWidth',1.5);
                                        end
                                    else
                                        if pV<pcritical(1) & isSignificant==1
                                            hline(yl(1)-0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNamesVerif{rmEffect(nRm)})(i,:),'lineWidth',1.5);
                                        else
                                            hline(yl(1)-0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNamesVerif{rmEffect(nRm)})(i,:),'lineWidth',1.5);
                                        end
                                    end
                                    nSignificant=nSignificant+1;
                                end
                            end
                        end


                        yl=yl*1.05.^max(nSignificant);
                        set(ax,'ylim',[min(yl) max(yl)])

                        print('-dtiff',['-r' num2str(imageResolution)], fullfile(saveDir, 'Lines', condNames{rmEffect(nRm)}, [condNames{indEffect(nInd)} ' = ' allMod.(condNamesVerif{indEffect(nInd)}){nMod}]))
                        close
                        clear ax yl dataMeans dataSD

                    end
                end
            end
        end

        % Text
        for nRm=1:numel(effectRM)
            for nInd=1:numel(indEffect)
                for nMod=1:numel(allMod.(condNamesVerif{indEffect(nInd)}))

                    f=figure('units','centimeters','position',[0 0 6+4*numel(modalitiesRM{nRm}) 4+9/16*4*numel(modalitiesRM{nRm})],'visible','off');

                    for x=1:numel(modalitiesRM{nRm})
                        dataMeans(x)=nanmean(nanmean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(:,col4means{nRm}(x,:)),2));
                        dataSD(x)=nanstd(nanmean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(:,col4means{nRm}(x,:)),2));
                    end

                    for x=1:numel(modalitiesRM{nRm})

                        h=bar(x,dataMeans(x)); hold on
                        h(1).FaceColor=colors{rmEffect(nRm)}(x,:);

                        if plotSD==1
                            if dataMeans(x)>=0
                                SD(1)=0;
                                SD(2)=dataSD(x);
                            else
                                SD(1)=dataSD(x);
                                SD(2)=0;
                            end
                            errorbar(x,dataMeans(x),SD(1),SD(2),'k','LineStyle','none')
                        end
                    end

                    if indivLines==1
                        for xl=1:numel(modalitiesRM{nRm})
                            selectIndivLines(:,xl)=mean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(:,col4means{nRm}(xl,:)),2);
                        end
                        xl=1:numel(modalitiesRM{nRm});
                        plot(xl,selectIndivLines,'k--'); hold on
                        scatter(xl,selectIndivLines,'k+')
                    end
                    clear selectIndivLines

                    for x=1:numel(modalitiesRM{nRm})
                        if abs(dataMeans(x))<1
                            text(x,0.5*(dataMeans(x)),sprintf('%.3f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                        elseif abs(dataMeans(x))<10
                            text(x,0.5*(dataMeans(x)),sprintf('%.2f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                        elseif abs(dataMeans(x))<100
                            text(x,0.5*(dataMeans(x)),sprintf('%.1f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                        else
                            text(x,0.5*(dataMeans(x)),sprintf('%.0f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                        end

                    end
                    xticks(1:numel(modalitiesRM{nRm}))
                    xticklabels(modalitiesRM{nRm})
                    xlabel(condNames{rmEffect(nRm)})
                    ylabel(units)
                    yl=ylim;
                    box off
                    ax=gca;
                    ax.XGrid='off';
                    ax.YGrid='on';

                    nAovInt=findcolExact(aov.Effect,[condNamesVerif{indEffect(nInd)} ':' condNamesVerif{rmEffect(nRm)}]);
                    nAov=findcolExact(aov.Effect, verifFieldName(condNames{rmEffect(nRm)}));
                    pINT=aov{nAovInt,6};
                    pMAIN=aov{nAov,6};

                    amp=[max(yl)-min(yl)];
                    isSignificant=0;

                    pValues=ones(1,size(postHoc.(condNamesVerif{rmEffect(nRm)}),1));
                    pSelected=0;
                    if pINT<pcritical(1)
                        isSignificant=1;
                        rows=findcolExact(postHoc.([condNamesVerif{rmEffect(nRm)} 'By' condNamesVerif{indEffect(nInd)}]){:,1},allMod.(condNamesVerif{indEffect(nInd)}){nMod});
                        pValues=postHoc.([condNamesVerif{rmEffect(nRm)} 'By' condNamesVerif{indEffect(nInd)}]){rows,6};
                        pSelected=1;
                    end
                    %                 if pMAIN<pcritical(1) &  pSelected==0
                    %                     isSignificant=1;
                    %                     pValues=postHoc.(condNamesVerif{rmEffect(nRm)}){:,5};
                    %                     pSelected=1;
                    %                 end
                    if pINT<pcritical(end) &  pSelected==0
                        rows=findcolExact(postHoc.([condNamesVerif{rmEffect(nRm)} 'By' condNamesVerif{indEffect(nInd)}]){:,1},allMod.(condNamesVerif{indEffect(nInd)}){nMod});
                        pValues=postHoc.([condNamesVerif{rmEffect(nRm)} 'By' condNamesVerif{indEffect(nInd)}]){rows,6};
                        pSelected=1;
                    end
                    %                 if pMAIN<pcritical(end) &  pSelected==0
                    %                     pValues=postHoc.(condNamesVerif{rmEffect(nRm)}){:,5};
                    %                     pSelected=1;
                    %                 end

                    if plotSD==1
                        dataMeans=dataMeans+sign(dataMeans).*dataSD;
                    end

                    nColSignificant=ones(1,numel(dataMeans));
                    nSignificant=1;
                    if pSelected==1
                        for i=1:numel(pValues)
                            pV=pValues(i);
                            if pV<=pcritical(end)
                                if abs(yl(2))>abs(yl(1))
                                    if abs(dataMeans(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1)))>abs(dataMeans(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2)))
                                        addPvalue(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1), yl(2)+0.055*nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1))*amp, pV, pcritical, colors{rmEffect(nRm)}(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2),:))
                                        nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1))=nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1))+1;
                                    else
                                        addPvalue(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2), yl(2)+0.055*nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2))*amp, pV, pcritical, colors{rmEffect(nRm)}(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1),:))
                                        nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2))=nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2))+1;

                                    end
                                else
                                    if abs(dataMeans(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1)))>abs(dataMeans(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2)))
                                        addPvalue(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1), yl(1)-0.055*nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1))*amp, pV, pcritical, colors{rmEffect(nRm)}(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2),:))
                                        nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1))=nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1))+1;
                                    else
                                        addPvalue(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2), yl(1)-0.055*nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2))*amp, pV, pcritical, colors{rmEffect(nRm)}(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1),:))
                                        nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2))=nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2))+1;
                                    end
                                end
                                nSignificant=nSignificant+1;
                            end
                        end
                    end

                    yl=yl*1.05.^max(nColSignificant);
                    set(ax,'ylim',[min(yl) max(yl)])

                    print('-dtiff',['-r' num2str(imageResolution)], fullfile(saveDir, 'Text', condNames{rmEffect(nRm)}, [condNames{indEffect(nInd)} ' = ' allMod.(condNamesVerif{indEffect(nInd)}){nMod}]))
                    close
                    clear ax yl dataMeans dataSD

                end
            end
        end
    end

    %% Main effect with double IND interaction
    if numel(rmEffect)>0 & numel(indEffect)>1
        % Lines
        if statsLines
            for nRm=1:numel(effectRM)
                for nInd1=1:numel(indEffect)
                    for nInd2=1:numel(indEffect)
                        if nInd2>nInd1

                            for nMod1=1:numel(allMod.(condNamesVerif{indEffect(nInd1)}))
                                for nMod2=1:numel(allMod.(condNamesVerif{indEffect(nInd2)}))

                                    f=figure('units','centimeters','position',[0 0 6+4*numel(modalitiesRM{nRm}) 4+9/16*4*numel(modalitiesRM{nRm})],'visible','off');

                                    for x=1:numel(modalitiesRM{nRm})
                                        dataMeans(x)=nanmean(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nMod1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2})(:,col4means{nRm}(x,:)),2));
                                        dataSD(x)=nanstd(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nMod1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2})(:,col4means{nRm}(x,:)),2));
                                    end

                                    for x=1:numel(modalitiesRM{nRm})

                                        h=bar(x,dataMeans(x)); hold on
                                        h(1).FaceColor=colors{rmEffect(nRm)}(x,:);

                                        if plotSD==1
                                            if dataMeans(x)>=0
                                                SD(1)=0;
                                                SD(2)=dataSD(x);
                                            else
                                                SD(1)=dataSD(x);
                                                SD(2)=0;
                                            end
                                            errorbar(x,dataMeans(x),SD(1),SD(2),'k','LineStyle','none')
                                        end

                                        if indivLines==1
                                            xl=1:numel(modalitiesRM{nRm});
                                            plot(xl,data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(:,col4means{nRm}),'k--'); hold on
                                            scatter(xl,data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(:,col4means{nRm}),'k+')
                                        end

                                        if abs(dataMeans(x))<1
                                            text(x,0.5*(dataMeans(x)),sprintf('%.3f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                        elseif abs(dataMeans(x))<10
                                            text(x,0.5*(dataMeans(x)),sprintf('%.2f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                        elseif abs(dataMeans(x))<100
                                            text(x,0.5*(dataMeans(x)),sprintf('%.1f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                        else
                                            text(x,0.5*(dataMeans(x)),sprintf('%.0f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                        end

                                    end
                                    xticks(1:numel(modalitiesRM{nRm}))
                                    xticklabels(modalitiesRM{nRm})
                                    xlabel(condNames{rmEffect(nRm)})
                                    ylabel(units)
                                    yl=ylim;
                                    box off
                                    ax=gca;
                                    ax.XGrid='off';
                                    ax.YGrid='on';

                                    nAovInt3=findcolExact(aov.Effect,[condNamesVerif{indEffect(nInd1)} ':' condNamesVerif{indEffect(nInd2)} ':' condNamesVerif{rmEffect(nRm)}]);
                                    if isempty(nAovInt3)
                                        nAovInt3=findComboStrings(aov.Effect,{condNamesVerif{indEffect(nInd1)}, condNamesVerif{indEffect(nInd2)}, condNamesVerif{rmEffect(nRm)}});
                                    end
                                    nAovInt2=findcolExact(aov.Effect,[condNamesVerif{rmEffect(nRm)} ':' condNamesVerif{indEffect(nInd1)}]);
                                    if isempty(nAovInt2)
                                        nAovInt2=findcolExact(aov.Effect,[condNamesVerif{rmEffect(nInd1)} ':' condNamesVerif{indEffect(nRm)}]);
                                    end
                                    nAovInt=findcolExact(aov.Effect,[condNamesVerif{rmEffect(nRm)} ':' condNamesVerif{indEffect(nInd2)}]);
                                    if isempty(nAovInt)
                                        nAovInt=findcolExact(aov.Effect,[condNamesVerif{rmEffect(nInd2)} ':' condNamesVerif{indEffect(nRm)}]);
                                    end
                                    nAov=findcolExact(aov.Effect, verifFieldName(condNames{rmEffect(nRm)}));
                                    pINT3=aov{nAovInt3,6};
                                    pINT2=aov{nAovInt2,6};
                                    pINT=aov{nAovInt,6};
                                    pMAIN=aov{nAov,6};

                                    amp=[max(yl)-min(yl)];
                                    isSignificant=0;

                                    pValues=ones(1,size(postHoc.(condNamesVerif{rmEffect(nRm)}),1));
                                    pSelected=0;
                                    if pINT3<pcritical(1)
                                        isSignificant=1;
                                        phCut1=findcolExact(postHoc.([condNamesVerif{rmEffect(nRm)} 'By' condNamesVerif{indEffect(nInd1)} 'By' condNamesVerif{indEffect(nInd2)}]){:,1},allMod.(condNamesVerif{indEffect(nInd1)}){nMod1});
                                        phCut2=findcolExact(postHoc.([condNamesVerif{rmEffect(nRm)} 'By' condNamesVerif{indEffect(nInd1)} 'By' condNamesVerif{indEffect(nInd2)}]){:,2},allMod.(condNamesVerif{indEffect(nInd2)}){nMod2});
                                        phCut=intersect(phCut1, phCut2);
                                        pValues=postHoc.([condNamesVerif{rmEffect(nRm)} 'By' condNamesVerif{indEffect(nInd1)} 'By' condNamesVerif{indEffect(nInd2)}]){phCut,7};
                                        pSelected=1;
                                    end
                                    %                                 if any([pINT<pcritical(1) pINT2<pcritical(1)]) &  pSelected==0
                                    %                                     isSignificant=1;
                                    %                                     pSelected=1;
                                    %                                     if pINT<pINT2
                                    %                                         phCut=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd2)}]){:,1}, allMod.(condNamesVerif{indEffect(nInd2)}){nMod2});
                                    %                                         pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd2)}]){phCut,6};
                                    %                                     else
                                    %                                         phCut=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)}]){:,1}, allMod.(condNamesVerif{indEffect(nInd1)}){nMod1});
                                    %                                         pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)}]){phCut,6};
                                    %                                     end
                                    %                                 end
                                    %                                 if pMAIN<pcritical(1) &  pSelected==0
                                    %                                     isSignificant=1;
                                    %                                     pValues=postHoc.(condNamesVerif{rmEffect(nRm)}){:,5};
                                    %                                     pSelected=1;
                                    %                                 end
                                    if pINT3<pcritical(end) &  pSelected==0
                                        phCut1=findcolExact(postHoc.([condNamesVerif{rmEffect(nRm)} 'By' condNamesVerif{indEffect(nInd1)} 'By' condNamesVerif{indEffect(nInd2)}]){:,1},allMod.(condNamesVerif{indEffect(nInd1)}){nMod1});
                                        phCut2=findcolExact(postHoc.([condNamesVerif{rmEffect(nRm)} 'By' condNamesVerif{indEffect(nInd1)} 'By' condNamesVerif{indEffect(nInd2)}]){:,2},allMod.(condNamesVerif{indEffect(nInd2)}){nMod2});
                                        phCut=intersect(phCut1, phCut2);
                                        pValues=postHoc.([condNamesVerif{rmEffect(nRm)} 'By' condNamesVerif{indEffect(nInd1)} 'By' condNamesVerif{indEffect(nInd2)}]){phCut,7};
                                        pSelected=1;
                                    end
                                    %                                 if any([pINT<pcritical(end) pINT2<pcritical(end)]) &  pSelected==0
                                    %                                     pSelected=1;
                                    %                                     if pINT<pINT2
                                    %                                         phCut=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd2)}]){:,1}, allMod.(condNamesVerif{indEffect(nInd2)}){nMod2});
                                    %                                         pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd2)}]){phCut,6};
                                    %                                     else
                                    %                                         phCut=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)}]){:,1}, allMod.(condNamesVerif{indEffect(nInd1)}){nMod1});
                                    %                                         pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)}]){phCut,6};
                                    %                                     end
                                    %                                 end
                                    %                                 if pMAIN<pcritical(end) &  pSelected==0
                                    %                                     pValues=postHoc.(condNamesVerif{rmEffect(nRm)}){:,5};
                                    %                                     pSelected=1;
                                    %                                 end

                                    if plotSD==1
                                        dataMeans=dataMeans+sign(dataMeans).*dataSD;
                                    end

                                    nSignificant=1;
                                    if pSelected==1
                                        nSignificant=1;
                                        for i=1:numel(pValues)
                                            pV=pValues(i);
                                            if pV<=pcritical(end)
                                                if abs(yl(2))>abs(yl(1))
                                                    if pV<pcritical(1) & isSignificant==1
                                                        hline(yl(2)+0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNamesVerif{rmEffect(nRm)})(i,:),'lineWidth',1.5);
                                                    else
                                                        hline(yl(2)+0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNamesVerif{rmEffect(nRm)})(i,:),'lineWidth',1.5);
                                                    end
                                                else
                                                    if pV<pcritical(1) & isSignificant==1
                                                        hline(yl(1)-0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNamesVerif{rmEffect(nRm)})(i,:),'lineWidth',1.5);
                                                    else
                                                        hline(yl(1)-0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNamesVerif{rmEffect(nRm)})(i,:),'lineWidth',1.5);
                                                    end
                                                end
                                                nSignificant=nSignificant+1;
                                            end
                                        end
                                    end


                                    yl=yl*1.05.^max(nSignificant);
                                    set(ax,'ylim',[min(yl) max(yl)])

                                    print('-dtiff',['-r' num2str(imageResolution)], fullfile(saveDir, 'Lines', condNames{rmEffect(nRm)}, [[condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}] ' = ' allMod.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){nMod}]))
                                    close
                                    clear ax yl dataMeans dataSD

                                end
                            end
                        end
                    end
                end
            end
        end

        % Text
        for nRm=1:numel(effectRM)
            for nInd1=1:numel(indEffect)
                for nInd2=1:numel(indEffect)
                    if nInd2>nInd1

                        for nMod1=1:numel(allMod.(condNamesVerif{indEffect(nInd1)}))
                            for nMod2=1:numel(allMod.(condNamesVerif{indEffect(nInd2)}))

                                f=figure('units','centimeters','position',[0 0 6+4*numel(modalitiesRM{nRm}) 4+9/16*4*numel(modalitiesRM{nRm})],'visible','off');

                                for x=1:numel(modalitiesRM{nRm})
                                    dataMeans(x)=nanmean(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nMod1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2})(:,col4means{nRm}(x,:)),2));
                                    dataSD(x)=nanstd(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nMod1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2})(:,col4means{nRm}(x,:)),2));
                                end

                                for x=1:numel(modalitiesRM{nRm})

                                    h=bar(x,dataMeans(x)); hold on
                                    h(1).FaceColor=colors{rmEffect(nRm)}(x,:);

                                    if plotSD==1
                                        if dataMeans(x)>=0
                                            SD(1)=0;
                                            SD(2)=dataSD(x);
                                        else
                                            SD(1)=dataSD(x);
                                            SD(2)=0;
                                        end
                                        errorbar(x,dataMeans(x),SD(1),SD(2),'k','LineStyle','none')
                                    end
                                end

                                if indivLines==1
                                    xl=1:numel(modalitiesRM{nRm});
                                    plot(xl,data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nMod1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2})(:,col4means{nRm}),'k--'); hold on
                                    scatter(xl,data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nMod1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2})(:,col4means{nRm}),'k+')
                                end

                                for x=1:numel(modalitiesRM{nRm})
                                    if abs(dataMeans(x))<1
                                        text(x,0.5*(dataMeans(x)),sprintf('%.3f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                    elseif abs(dataMeans(x))<10
                                        text(x,0.5*(dataMeans(x)),sprintf('%.2f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                    elseif abs(dataMeans(x))<100
                                        text(x,0.5*(dataMeans(x)),sprintf('%.1f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                    else
                                        text(x,0.5*(dataMeans(x)),sprintf('%.0f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                    end

                                end
                                xticks(1:numel(modalitiesRM{nRm}))
                                xticklabels(modalitiesRM{nRm})
                                xlabel(condNames{rmEffect(nRm)})
                                ylabel(units)
                                yl=ylim;
                                box off
                                ax=gca;
                                ax.XGrid='off';
                                ax.YGrid='on';

                                nAovInt3=findcolExact(aov.Effect,[condNames{indEffect(nInd1)} ':' condNames{indEffect(nInd2)} ':' condNames{rmEffect(nRm)}]);
                                if isempty(nAovInt3)
                                    nAovInt3=findComboStrings(aov.Effect,{condNames{indEffect(nInd1)}, condNames{indEffect(nInd2)}, condNames{rmEffect(nRm)}});
                                end
                                %                             nAovInt2=findcolExact(aov.Effect,[condNames{rmEffect(nRm)} ':' condNames{indEffect(nInd1)}]);
                                %                             if isempty(nAovInt2)
                                %                                 nAovInt2=findcolExact(aov.Effect,[condNames{rmEffect(nInd1)} ':' condNames{indEffect(nRm)}]);
                                %                             end
                                %                             nAovInt=findcolExact(aov.Effect,[condNames{rmEffect(nRm)} ':' condNames{indEffect(nInd2)}]);
                                %                             if isempty(nAovInt)
                                %                                 nAovInt=findcolExact(aov.Effect,[condNames{rmEffect(nInd2)} ':' condNames{indEffect(nRm)}]);
                                %                             end
                                %                             nAov=findcolExact(aov.Effect, verifFieldName(condNames{rmEffect(nRm)});
                                pINT3=aov{nAovInt3,6};
                                %                             pINT2=aov{nAovInt2,6};
                                %                             pINT=aov{nAovInt,6};
                                %                             pMAIN=aov{nAov,6};

                                amp=[max(yl)-min(yl)];
                                isSignificant=0;

                                pValues=ones(1,size(postHoc.(condNamesVerif{rmEffect(nRm)}),1));
                                pSelected=0;
                                if pINT3<pcritical(1)
                                    isSignificant=1;
                                    phCut1=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,1},allMod.(condNamesVerif{indEffect(nInd1)}){nMod1});
                                    phCut2=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,2},allMod.(condNamesVerif{indEffect(nInd2)}){nMod2});
                                    phCut=intersect(phCut1, phCut2);
                                    pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){phCut,7};
                                    pSelected=1;
                                end
                                %                             if any([pINT<pcritical(1) pINT2<pcritical(1)]) &  pSelected==0
                                %                                 isSignificant=1;
                                %                                 pSelected=1;
                                %                                 if pINT<pINT2
                                %                                     phCut=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd2)}]){:,1}, allMod.(condNamesVerif{indEffect(nInd2)}){nMod2});
                                %                                     pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd2)}]){phCut,6};
                                %                                 else
                                %                                     phCut=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)}]){:,1}, allMod.(condNamesVerif{indEffect(nInd1)}){nMod1});
                                %                                     pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)}]){phCut,6};
                                %                                 end
                                %                             end
                                %                             if pMAIN<pcritical(1) &  pSelected==0
                                %                                 isSignificant=1;
                                %                                 pValues=postHoc.(condNamesVerif{rmEffect(nRm)}){:,5};
                                %                                 pSelected=1;
                                %                             end
                                if pINT3<pcritical(end) &  pSelected==0
                                    phCut1=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,1},allMod.(condNamesVerif{indEffect(nInd1)}){nMod1});
                                    phCut2=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,2},allMod.(condNamesVerif{indEffect(nInd2)}){nMod2});
                                    phCut=intersect(phCut1, phCut2);
                                    pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){phCut,7};
                                    pSelected=1;
                                end
                                %                             if any([pINT<pcritical(end) pINT2<pcritical(end)]) &  pSelected==0
                                %                                 pSelected=1;
                                %                                 if pINT<pINT2
                                %                                     phCut=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd2)}]){:,1}, allMod.(condNamesVerif{indEffect(nInd2)}){nMod2});
                                %                                     pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd2)}]){phCut,6};
                                %                                 else
                                %                                     phCut=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)}]){:,1}, allMod.(condNamesVerif{indEffect(nInd1)}){nMod1});
                                %                                     pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)}]){phCut,6};
                                %                                 end
                                %                             end
                                %                             if pMAIN<pcritical(end) &  pSelected==0
                                %                                 pValues=postHoc.(condNamesVerif{rmEffect(nRm)}){:,5};
                                %                                 pSelected=1;
                                %                             end

                                if plotSD==1
                                    dataMeans=dataMeans+sign(dataMeans).*dataSD;
                                end

                                nColSignificant=ones(1,numel(dataMeans));
                                nSignificant=1;
                                if pSelected==1
                                    for i=1:numel(pValues)
                                        pV=pValues(i);
                                        if pV<=pcritical(end)
                                            if abs(yl(2))>abs(yl(1))
                                                if abs(dataMeans(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1)))>abs(dataMeans(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2)))
                                                    addPvalue(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1), yl(2)+0.055*nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1))*amp, pV, pcritical, colors{rmEffect(nRm)}(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2),:))
                                                    nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1))=nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1))+1;
                                                else
                                                    addPvalue(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2), yl(2)+0.055*nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2))*amp, pV, pcritical, colors{rmEffect(nRm)}(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1),:))
                                                    nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2))=nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2))+1;

                                                end
                                            else
                                                if abs(dataMeans(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1)))>abs(dataMeans(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2)))
                                                    addPvalue(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1), yl(1)-0.055*nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1))*amp, pV, pcritical, colors{rmEffect(nRm)}(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2),:))
                                                    nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1))=nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1))+1;
                                                else
                                                    addPvalue(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2), yl(1)-0.055*nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2))*amp, pV, pcritical, colors{rmEffect(nRm)}(order4ES.(condNamesVerif{rmEffect(nRm)})(i,1),:))
                                                    nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2))=nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm)})(i,2))+1;
                                                end
                                            end
                                            nSignificant=nSignificant+1;
                                        end
                                    end
                                end

                                yl=yl*1.05.^max(nColSignificant);
                                set(ax,'ylim',[min(yl) max(yl)])

                                print('-dtiff',['-r' num2str(imageResolution)], fullfile(saveDir, 'Text', condNames{rmEffect(nRm)}, [condNames{indEffect(nInd1)} ' By ' condNames{indEffect(nInd2)} ' = ' allMod.(condNamesVerif{indEffect(nInd1)}){nMod1} ' ' allMod.(condNamesVerif{indEffect(nInd2)}){nMod2}]))
                                close
                                clear ax yl dataMeans dataSD

                            end
                        end
                    end
                end
            end
        end
    end

    %% Interaction of RM effect with simple IND effect
    if numel(rmEffect)>1 & numel(indEffect)>0
        % Lines
        if statsLines
            for nRm1=1:numel(modalitiesRM)
                for nRm2=1:numel(modalitiesRM)
                    if nRm1~=nRm2
                        for nInd=1:numel(indEffect)
                            for nMod=1:numel(allMod.(condNamesVerif{indEffect(nInd)}))

                                f=figure('units','centimeters','position',[0 0 6+4*numel(allModalities{rmEffect(nRm1)}) 4+numel(cond4effect{rmEffect(nRm2)})*numel(cond4effect{rmEffect(nRm1)})*9/16*4],'visible','off');

                                for nModRm=1:numel(cond4effect{rmEffect(nRm2)})
                                    dataMeans=nanmean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(:,col4means{nRm1}(:,nModRm)));
                                    dataSD=nanstd(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(:,col4means{nRm1}(:,nModRm)));

                                    subplot(numel(cond4effect{rmEffect(nRm2)}),1,nModRm);
                                    for x=1:numel(cond4effect{rmEffect(nRm1)})
                                        h=bar(x,dataMeans(x)); hold on
                                        h(1).FaceColor=colors{rmEffect(nRm1)}(x,:);

                                        if plotSD==1
                                            if dataMeans(x)>=0
                                                SD(1)=0;
                                                SD(2)=dataSD(x);
                                            else
                                                SD(1)=dataSD(x);
                                                SD(2)=0;
                                            end
                                            errorbar(x,dataMeans(x),SD(1),SD(2),'k','LineStyle','none')
                                        end
                                    end

                                    if indivLines==1
                                        xl=1:numel(modalitiesRM{nRm1});
                                        plot(xl,data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(:,col4means{nRm1}(:,nModRm)),'k--'); hold on
                                        scatter(xl,data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(:,col4means{nRm1}(:,nModRm)),'k+')
                                    end

                                    for x=1:numel(cond4effect{rmEffect(nRm1)})
                                        if abs(dataMeans(x))<1
                                            text(x,0.5*(dataMeans(x)),sprintf('%.3f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                        elseif abs(dataMeans(x))<10
                                            text(x,0.5*(dataMeans(x)),sprintf('%.2f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                        elseif abs(dataMeans(x))<100
                                            text(x,0.5*(dataMeans(x)),sprintf('%.1f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                        else
                                            text(x,0.5*(dataMeans(x)),sprintf('%.0f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                        end
                                    end

                                    xticks(1:numel(cond4effect{rmEffect(nRm1)}))
                                    xticklabels(cond4effect{rmEffect(nRm1)})
                                    if nModRm==numel(cond4effect{rmEffect(nRm2)})
                                        xlabel(condNames{rmEffect(nRm1)})
                                    end
                                    ylabel(units)
                                    box off
                                    ax{nModRm}=gca;
                                    ax{nModRm}.XGrid='off';
                                    ax{nModRm}.YGrid='on';
                                    yl(:,nModRm)=ylim;
                                    title(cond4effect{rmEffect(nRm2)}{nModRm})

                                end

                                nAovInt3=findcolExact(aov.Effect,[condNames{indEffect(nInd)} ':' condNames{rmEffect(nRm1)} ':' condNames{rmEffect(nRm2)}]);
                                if isempty(nAovInt3)
                                    nAovInt3=findComboStrings(aov.Effect,{condNames{indEffect(nInd)}, condNames{rmEffect(nRm1)}, condNames{rmEffect(nRm2)}});
                                end
                                nAovInt2=findcolExact(aov.Effect,[condNames{indEffect(nInd)} ':' condNames{rmEffect(nRm1)}]);
                                if isempty(nAovInt2)
                                    nAovInt2=findcolExact(aov.Effect,[condNames{rmEffect(nRm1)} ':' condNames{indEffect(nInd)}]);
                                end
                                nAovInt=findcolExact(aov.Effect,[condNames{rmEffect(nRm1)} ':' condNames{rmEffect(nRm2)}]);
                                if isempty(nAovInt)
                                    nAovInt=findcolExact(aov.Effect,[condNames{rmEffect(nRm2)} ':' condNames{rmEffect(nRm1)}]);
                                end

                                nAov=findcolExact(aov.Effect, verifFieldName(condNames{rmEffect(nRm1)}));
                                pINT3=aov{nAovInt3,6};
                                pINT2=aov{nAovInt2,6};
                                pINT=aov{nAovInt,6};
                                pMAIN=aov{nAov,6};

                                amp=[max(max(yl))-min(min(yl))];
                                for nModRm=1:numel(cond4effect{rmEffect(nRm2)})
                                    isSignificant=0;
                                    set(f,'CurrentAxes',ax{nModRm});
                                    pValues=ones(1,size(postHoc.(condNamesVerif{rmEffect(nRm1)}),1));
                                    pSelected=0;

                                    if pINT3<pcritical(1)
                                        isSignificant=1;
                                        phCut1=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd)}]){:,1}, allMod.(condNamesVerif{indEffect(nInd)}){nMod});
                                        phCut2=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd)}]){:,2}, allMod.(condNamesVerif{rmEffect(nRm2)}){nModRm});
                                        phCut=intersect(phCut1, phCut2);
                                        pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd)}]){phCut,7};
                                        pSelected=1;
                                    end
                                    %                                 if any([pINT<pcritical(1) pINT2<pcritical(1)]) &  pSelected==0
                                    %                                     isSignificant=1;
                                    %                                     pSelected=1;
                                    %                                     if pINT<pINT2
                                    %                                         phCut=findcolExact(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}]){:,1}, allMod.(condNamesVerif{rmEffect(nRm2)}){nModRm});
                                    %                                         pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}]){phCut,6};
                                    %                                     else
                                    %                                         phCut=findcolExact(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{indEffect(nInd)}]){:,1}, allMod.(condNamesVerif{indEffect(nInd)}){nMod});
                                    %                                         pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{indEffect(nInd)}]){phCut,6};
                                    %                                     end
                                    %                                 end
                                    %                                 if pMAIN<pcritical(1) &  pSelected==0
                                    %                                     isSignificant=1;
                                    %                                     pValues=postHoc.(condNamesVerif{rmEffect(nRm1)}){:,5};
                                    %                                     pSelected=1;
                                    %                                 end
                                    if pINT3<pcritical(end) &  pSelected==0
                                        phCut1=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd)}]){:,1}, allMod.(condNamesVerif{indEffect(nInd)}){nMod});
                                        phCut2=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd)}]){:,2}, allMod.(condNamesVerif{rmEffect(nRm2)}){nModRm});
                                        phCut=intersect(phCut1, phCut2);
                                        pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd)}]){phCut,7};
                                        pSelected=1;
                                    end
                                    %                                 if any([pINT<pcritical(end) pINT2<pcritical(end)]) &  pSelected==0
                                    %                                     pSelected=1;
                                    %                                     if pINT<pINT2
                                    %                                         phCut=findcolExact(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}]){:,1}, allMod.(condNamesVerif{rmEffect(nRm2)}){nModRm});
                                    %                                         pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}]){phCut,6};
                                    %                                     else
                                    %                                         phCut=findcolExact(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{indEffect(nInd)}]){:,1}, allMod.(condNamesVerif{indEffect(nInd)}){nMod});
                                    %                                         pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{indEffect(nInd)}]){phCut,6};
                                    %                                     end
                                    %                                 end
                                    %                                 if pMAIN<pcritical(end) &  pSelected==0
                                    %                                     pValues=postHoc.(condNamesVerif{rmEffect(nRm1)}){:,5};
                                    %                                     pSelected=1;
                                    %                                 end

                                    if pSelected==1
                                        nSignificant=1;
                                        for i=1:numel(pValues)
                                            pV=pValues(i);
                                            if pV<=pcritical(end)
                                                if yl(2,nModRm)>0
                                                    if pV<pcritical(1) & isSignificant==1
                                                        hline(yl(2,nModRm)+0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNamesVerif{rmEffect(nRm1)})(i,:),'lineWidth',1.5);
                                                    else
                                                        hline(yl(2,nModRm)+0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNamesVerif{rmEffect(nRm1)})(i,:),'lineWidth',1.5);
                                                    end
                                                else
                                                    if pV<pcritical(1) & isSignificant==1
                                                        hline(yl(1,nModRm)-0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNamesVerif{rmEffect(nRm1)})(i,:),'lineWidth',1.5);
                                                    else
                                                        hline(yl(1,nModRm)-0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNamesVerif{rmEffect(nRm1)})(i,:),'lineWidth',1.5);
                                                    end
                                                end
                                                nSignificant=nSignificant+1;
                                            end
                                        end
                                    end

                                    yl(:,nModRm)=ylim;

                                end

                                for nModRm=1:numel(ax)
                                    set(ax{nModRm},'ylim',1.05*[min(min(yl)) max(max(yl))])
                                end

                                print('-dtiff',['-r' num2str(imageResolution)],fullfile(saveDir, 'Lines', condNames{rmEffect(nRm1)}, [condNames{indEffect(nInd)} ' = ' allMod.(condNamesVerif{indEffect(nInd)}){nMod} ' By ' condNames{rmEffect(nRm2)}]))
                                close
                                clear ax yl dataMeans dataSD

                            end
                        end
                    end
                end
            end
        end

        % text
        for nRm1=1:numel(modalitiesRM)
            for nRm2=1:numel(modalitiesRM)
                if nRm1~=nRm2
                    for nInd=1:numel(indEffect)
                        for nMod=1:numel(allMod.(condNamesVerif{indEffect(nInd)}))

                            f=figure('units','centimeters','position',[0 0 6+4*numel(allModalities{rmEffect(nRm1)}) 4+numel(cond4effect{rmEffect(nRm2)})*numel(cond4effect{rmEffect(nRm1)})*9/16*4],'visible','off');

                            for nModRm=1:numel(cond4effect{rmEffect(nRm2)})
                                dataMeans=nanmean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(:,col4means{nRm1}(:,nModRm)));
                                dataSD=nanstd(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(:,col4means{nRm1}(:,nModRm)));

                                subplot(numel(cond4effect{rmEffect(nRm2)}),1,nModRm);

                                for x=1:numel(cond4effect{rmEffect(nRm1)})
                                    h=bar(x,dataMeans(x)); hold on
                                    h(1).FaceColor=colors{rmEffect(nRm1)}(x,:);

                                    if plotSD==1
                                        if dataMeans(x)>=0
                                            SD(1)=0;
                                            SD(2)=dataSD(x);
                                        else
                                            SD(1)=dataSD(x);
                                            SD(2)=0;
                                        end
                                        errorbar(x,dataMeans(x),SD(1),SD(2),'k','LineStyle','none')
                                    end
                                end

                                if indivLines==1
                                    xl=1:numel(modalitiesRM{nRm1});
                                    plot(xl,data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(:,col4means{nRm1}(:,nModRm)),'k--'); hold on
                                    scatter(xl,data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(:,col4means{nRm1}(:,nModRm)),'k+')
                                end

                                for x=1:numel(cond4effect{rmEffect(nRm1)})
                                    if abs(dataMeans(x))<1
                                        text(x,0.5*(dataMeans(x)),sprintf('%.3f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                    elseif abs(dataMeans(x))<10
                                        text(x,0.5*(dataMeans(x)),sprintf('%.2f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                    elseif abs(dataMeans(x))<100
                                        text(x,0.5*(dataMeans(x)),sprintf('%.1f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                    else
                                        text(x,0.5*(dataMeans(x)),sprintf('%.0f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                    end
                                end

                                xticks(1:numel(cond4effect{rmEffect(nRm1)}))
                                xticklabels(cond4effect{rmEffect(nRm1)})
                                if nModRm==numel(cond4effect{rmEffect(nRm2)})
                                    xlabel(condNames{rmEffect(nRm1)})
                                end
                                ylabel(units)
                                box off
                                ax{nModRm}=gca;
                                ax{nModRm}.XGrid='off';
                                ax{nModRm}.YGrid='on';
                                yl(:,nModRm)=ylim;
                                title(cond4effect{rmEffect(nRm2)}{nModRm})

                            end

                            nAovInt3=findcolExact(aov.Effect,[condNames{indEffect(nInd)} ':' condNames{rmEffect(nRm1)} ':' condNames{rmEffect(nRm2)}]);
                            if isempty(nAovInt3)
                                nAovInt3=findComboStrings(aov.Effect,{condNames{indEffect(nInd)}, condNames{rmEffect(nRm1)}, condNames{rmEffect(nRm2)}});
                            end
                            nAovInt2=findcolExact(aov.Effect,[condNames{indEffect(nInd)} ':' condNames{rmEffect(nRm1)}]);
                            if isempty(nAovInt2)
                                nAovInt2=findcolExact(aov.Effect,[condNames{rmEffect(nRm1)} ':' condNames{indEffect(nInd)}]);
                            end
                            nAovInt=findcolExact(aov.Effect,[condNames{rmEffect(nRm1)} ':' condNames{rmEffect(nRm2)}]);
                            if isempty(nAovInt)
                                nAovInt=findcolExact(aov.Effect,[condNames{rmEffect(nRm2)} ':' condNames{rmEffect(nRm1)}]);
                            end
                            nAov=findcolExact(aov.Effect, verifFieldName(condNames{rmEffect(nRm1)}));
                            pINT3=aov{nAovInt3,6};
                            pINT2=aov{nAovInt2,6};
                            pINT=aov{nAovInt,6};
                            pMAIN=aov{nAov,6};

                            amp=[max(max(yl))-min(min(yl))];
                            for nModRm=1:numel(cond4effect{rmEffect(nRm2)})

                                dataMeans=nanmean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(:,col4means{nRm1}(:,nModRm)));
                                dataSD=nanstd(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(:,col4means{nRm1}(:,nModRm)));

                                dataMeans4pv=dataMeans;
                                if plotSD==1
                                    dataMeans4pv=dataMeans+sign(dataMeans).*dataSD;
                                end
                                if indivLines==1
                                    dataMeans4pv=max(abs(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(:,col4means{nRm1}(:,nModRm)))).*sign(dataMeans);
                                end
                                if plotSD==1 && indivLines==1
                                    dataMeans4pv=sign(dataMeans).*max([abs(dataMeans+sign(dataMeans).*dataSD); max(abs(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(:,col4means{nRm1}(:,nModRm)))).*sign(dataMeans)]);
                                end


                                isSignificant=0;
                                set(f,'CurrentAxes',ax{nModRm});
                                pValues=ones(1,size(postHoc.(condNamesVerif{rmEffect(nRm1)}),1));
                                pSelected=0;
                                if pINT3<pcritical(1)
                                    isSignificant=1;
                                    phCut1=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd)}]){:,1}, allMod.(condNamesVerif{indEffect(nInd)}){nMod});
                                    phCut2=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd)}]){:,2}, allMod.(condNamesVerif{rmEffect(nRm2)}){nModRm});
                                    phCut=intersect(phCut1, phCut2);
                                    pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd)}]){phCut,7};
                                    pSelected=1;
                                end
                                %                             if any([pINT<pcritical(1) pINT2<pcritical(1)]) &  pSelected==0
                                %                                 isSignificant=1;
                                %                                 pSelected=1;
                                %                                 if pINT<pINT2
                                %                                     phCut=findcolExact(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}]){:,1}, allMod.(condNamesVerif{rmEffect(nRm2)}){nModRm});
                                %                                     pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}]){phCut,6};
                                %                                 else
                                %                                     phCut=findcolExact(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{indEffect(nInd)}]){:,1}, allMod.(condNamesVerif{indEffect(nInd)}){nMod});
                                %                                     pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{indEffect(nInd)}]){phCut,6};
                                %                                 end
                                %                             end
                                %                             if pMAIN<pcritical(1) &  pSelected==0
                                %                                 isSignificant=1;
                                %                                 pValues=postHoc.(condNamesVerif{rmEffect(nRm1)}){:,5};
                                %                                 pSelected=1;
                                %                             end
                                if pINT3<pcritical(end) &  pSelected==0
                                    phCut1=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd)}]){:,1}, allMod.(condNamesVerif{indEffect(nInd)}){nMod});
                                    phCut2=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd)}]){:,2}, allMod.(condNamesVerif{rmEffect(nRm2)}){nModRm});
                                    phCut=intersect(phCut1, phCut2);
                                    pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd)}]){phCut,7};
                                    pSelected=1;
                                end
                                %                             if any([pINT<pcritical(end) pINT2<pcritical(end)]) &  pSelected==0
                                %                                 pSelected=1;
                                %                                 if pINT<pINT2
                                %                                     phCut=findcolExact(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}]){:,1}, allMod.(condNamesVerif{rmEffect(nRm2)}){nModRm});
                                %                                     pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}]){phCut,6};
                                %                                 else
                                %                                     phCut=findcolExact(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{indEffect(nInd)}]){:,1}, allMod.(condNamesVerif{indEffect(nInd)}){nMod});
                                %                                     pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{indEffect(nInd)}]){phCut,6};
                                %                                 end
                                %                             end
                                %                             if pMAIN<pcritical(end) &  pSelected==0
                                %                                 pValues=postHoc.(condNamesVerif{rmEffect(nRm1)}){:,5};
                                %                                 pSelected=1;
                                %                             end

                                nSignificant=1;
                                nColSignificant=ones(1,numel(dataMeans));
                                if pSelected==1
                                    for i=1:numel(pValues)
                                        pV=pValues(i);
                                        if pV<=pcritical(end)
                                            if abs(yl(2))>abs(yl(1))
                                                if abs(dataMeans(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1)))>abs(dataMeans(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2)))
                                                    addPvalue(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1), dataMeans4pv(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1))+0.0225*numel(cond4effect{rmEffect(nRm2)})*nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1))*amp, pV, pcritical, colors{rmEffect(nRm1)}(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2),:))
                                                    nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1))=nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1))+1;
                                                else
                                                    addPvalue(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2), dataMeans(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2))+0.0225*numel(cond4effect{rmEffect(nRm2)})*nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2))*amp, pV, pcritical, colors{rmEffect(nRm1)}(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1),:))
                                                    nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2))=nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2))+1;
                                                end
                                            else
                                                if abs(dataMeans(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1)))>abs(dataMeans(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2)))
                                                    addPvalue(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1), dataMeans4pv(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1))-0.0225*numel(cond4effect{rmEffect(nRm2)})*nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1))*amp, pV, pcritical, colors{rmEffect(nRm1)}(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2),:))
                                                    nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1))=nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1))+1;
                                                else
                                                    addPvalue(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2), dataMeans4pv(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2))-0.0225*numel(cond4effect{rmEffect(nRm2)})*nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2))*amp, pV, pcritical, colors{rmEffect(nRm1)}(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1),:))
                                                    nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2))=nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2))+1;
                                                end
                                            end
                                            nSignificant=nSignificant+1;
                                        end
                                    end
                                end

                                yl(:,nModRm)=ylim*1.05.^max(nColSignificant);

                            end

                            for nModRm=1:numel(ax)
                                set(ax{nModRm},'ylim',[min(min(yl)) max(max(yl))])
                            end

                            print('-dtiff',['-r' num2str(imageResolution)],fullfile(saveDir, 'Text', condNames{rmEffect(nRm1)}, [condNames{indEffect(nInd)} ' = ' allMod.(condNamesVerif{indEffect(nInd)}){nMod} ' By ' condNames{rmEffect(nRm2)}]))
                            close
                            clear ax yl dataMeans dataSD

                        end
                    end
                end
            end
        end
    end

    %% Interaction of 2RM with 2IND
    if numel(rmEffect)>1 & numel(indEffect)>1
        % Lines
        if statsLines
            for nRm1=1:numel(modalitiesRM)
                for nRm2=1:numel(modalitiesRM)
                    if nRm1~=nRm2
                        for nInd1=1:numel(indEffect)
                            for nInd2=1:numel(indEffect)
                                if nInd2>nInd1
                                    for nModInd1=1:numel(allMod.(condNamesVerif{indEffect(nInd1)}))
                                        for nModInd2=1:numel(allMod.(condNamesVerif{indEffect(nInd2)}))

                                            f=figure('units','centimeters','position',[0 0 6+4*numel(allModalities{rmEffect(nRm1)}) 4+numel(cond4effect{rmEffect(nRm2)})*numel(cond4effect{rmEffect(nRm1)})*9/16*4],'visible','off');

                                            for nModRm=1:numel(cond4effect{rmEffect(nRm2)})
                                                dataMeans=nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nModInd1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nModInd2})(:,col4means{nRm1}(:,nModRm)));
                                                dataSD=nanstd(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nModInd1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nModInd2})(:,col4means{nRm1}(:,nModRm)));

                                                subplot(numel(cond4effect{rmEffect(nRm2)}),1,nModRm);
                                                for x=1:numel(cond4effect{rmEffect(nRm1)})
                                                    h=bar(x,dataMeans(x)); hold on
                                                    h(1).FaceColor=colors{rmEffect(nRm1)}(x,:);

                                                    if plotSD==1
                                                        if dataMeans(x)>=0
                                                            SD(1)=0;
                                                            SD(2)=dataSD(x);
                                                        else
                                                            SD(1)=dataSD(x);
                                                            SD(2)=0;
                                                        end
                                                        errorbar(x,dataMeans(x),SD(1),SD(2),'k','LineStyle','none')
                                                    end
                                                end

                                                if indivLines==1
                                                    xl=1:numel(modalitiesRM{nRm});
                                                    plot(xl,data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nModInd1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nModInd2})(:,col4means{nRm1}(:,nModRm)),'k--'); hold on
                                                    scatter(xl,data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nModInd1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nModInd2})(:,col4means{nRm1}(:,nModRm)),'k+')
                                                end

                                                for x=1:numel(cond4effect{rmEffect(nRm1)})
                                                    if abs(dataMeans(x))<1
                                                        text(x,0.5*(dataMeans(x)),sprintf('%.3f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                                    elseif abs(dataMeans(x))<10
                                                        text(x,0.5*(dataMeans(x)),sprintf('%.2f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                                    elseif abs(dataMeans(x))<100
                                                        text(x,0.5*(dataMeans(x)),sprintf('%.1f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                                    else
                                                        text(x,0.5*(dataMeans(x)),sprintf('%.0f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                                    end
                                                end

                                                xticks(1:numel(cond4effect{rmEffect(nRm1)}))
                                                xticklabels(cond4effect{rmEffect(nRm1)})
                                                if nModRm==numel(cond4effect{rmEffect(nRm2)})
                                                    xlabel(condNames{rmEffect(nRm1)})
                                                end
                                                ylabel(units)
                                                box off
                                                ax{nModRm}=gca;
                                                ax{nModRm}.XGrid='off';
                                                ax{nModRm}.YGrid='on';
                                                yl(:,nModRm)=ylim;
                                                title(cond4effect{rmEffect(nRm2)}{nModRm})

                                            end

                                            nAovInt4=findComboStrings(aov.Effect,{condNames{indEffect(nInd1)}, condNames{indEffect(nInd2)}, condNames{rmEffect(nRm1)}, condNames{rmEffect(nRm2)}});
                                            pINT4=aov{nAovInt4,6};

                                            amp=[max(max(yl))-min(min(yl))];

                                            for nModRm=1:numel(cond4effect{rmEffect(nRm2)})

                                                dataMeans=nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nModInd1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nModInd2})(:,col4means{nRm1}(:,nModRm)));
                                                dataSD=nanstd(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nModInd1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nModInd2})(:,col4means{nRm1}(:,nModRm)));

                                                if plotSD==1
                                                    dataMeans=dataMeans+sign(dataMeans).*dataSD;
                                                end

                                                isSignificant=0;
                                                set(f,'CurrentAxes',ax{nModRm});
                                                pValues=ones(1,size(postHoc.(condNamesVerif{rmEffect(nRm1)}),1));
                                                pSelected=0;
                                                if pINT4<pcritical(1)
                                                    isSignificant=1;
                                                    ph4=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]);
                                                    phCut1=findcol(ph4{:,1}, allMod.(condNamesVerif{indEffect(nInd1)}){nModInd1});
                                                    phCut2=findcol(ph4{:,2}, allMod.(condNamesVerif{indEffect(nInd2)}){nModInd2});
                                                    phCut3=findcol(ph4{:,3}, allMod.(condNamesVerif{rmEffect(nRm2)}){nModRm});
                                                    phCut4=intersect(phCut1, phCut2);
                                                    phCut=intersect(phCut4, phCut3);
                                                    pValues=ph4{phCut,8};
                                                    pSelected=1;
                                                end
                                                if pINT4<pcritical(end) & pSelected==0
                                                    ph4=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]);
                                                    phCut1=findcol(ph4{:,1}, allMod.(condNamesVerif{indEffect(nInd1)}){nModInd1});
                                                    phCut2=findcol(ph4{:,2}, allMod.(condNamesVerif{indEffect(nInd2)}){nModInd2});
                                                    phCut3=findcol(ph4{:,3}, allMod.(condNamesVerif{rmEffect(nRm2)}){nModRm});
                                                    phCut4=intersect(phCut1, phCut2);
                                                    phCut=intersect(phCut4, phCut3);
                                                    pValues=ph4{phCut,8};
                                                    pSelected=1;
                                                end
                                                if pSelected==1
                                                    nSignificant=1;
                                                    for i=1:numel(pValues)
                                                        pV=pValues(i);
                                                        if pV<=pcritical(end)
                                                            if yl(2,nModRm)>0
                                                                if pV<pcritical(1) & isSignificant==1
                                                                    hline(yl(2,nModRm)+0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNamesVerif{rmEffect(nRm1)})(i,:),'lineWidth',1.5);
                                                                else
                                                                    hline(yl(2,nModRm)+0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNamesVerif{rmEffect(nRm1)})(i,:),'lineWidth',1.5);
                                                                end
                                                            else
                                                                if pV<pcritical(1) & isSignificant==1
                                                                    hline(yl(1,nModRm)-0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNamesVerif{rmEffect(nRm1)})(i,:),'lineWidth',1.5);
                                                                else
                                                                    hline(yl(1,nModRm)-0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNamesVerif{rmEffect(nRm1)})(i,:),'lineWidth',1.5);
                                                                end
                                                            end
                                                            nSignificant=nSignificant+1;
                                                        end
                                                    end
                                                end

                                                yl(:,nModRm)=ylim;

                                            end

                                            for nModRm=1:numel(ax)
                                                set(ax{nModRm},'ylim',1.05*[min(min(yl)) max(max(yl))])
                                            end

                                            print('-dtiff',['-r' num2str(imageResolution)],fullfile(saveDir, 'Lines', condNames{rmEffect(nRm1)}, [condNames{indEffect(nInd1)} ' By ' condNames{indEffect(nInd2)} ' = ' allMod.(condNamesVerif{indEffect(nInd1)}){nModInd1} ' ' allMod.(condNamesVerif{indEffect(nInd2)}){nModInd2} ' By ' condNames{rmEffect(nRm2)}]))
                                            close
                                            clear ax yl dataMeans dataSD

                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        % text
        for nRm1=1:numel(modalitiesRM)
            for nRm2=1:numel(modalitiesRM)
                if nRm1~=nRm2
                    for nInd1=1:numel(indEffect)
                        for nInd2=1:numel(indEffect)
                            if nInd2>nInd1
                                for nModInd1=1:numel(allMod.(condNamesVerif{indEffect(nInd1)}))
                                    for nModInd2=1:numel(allMod.(condNamesVerif{indEffect(nInd2)}))

                                        f=figure('units','centimeters','position',[0 0 6+4*numel(allModalities{rmEffect(nRm1)}) 4+numel(cond4effect{rmEffect(nRm2)})*numel(cond4effect{rmEffect(nRm1)})*9/16*4],'visible','off');

                                        for nModRm=1:numel(cond4effect{rmEffect(nRm2)})
                                            dataMeans=nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nModInd1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nModInd2})(:,col4means{nRm1}(:,nModRm)));
                                            dataSD=nanstd(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nModInd1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nModInd2})(:,col4means{nRm1}(:,nModRm)));

                                            subplot(numel(cond4effect{rmEffect(nRm2)}),1,nModRm);
                                            for x=1:numel(cond4effect{rmEffect(nRm1)})
                                                h=bar(x,dataMeans(x)); hold on
                                                h(1).FaceColor=colors{rmEffect(nRm1)}(x,:);

                                                if plotSD==1
                                                    if dataMeans(x)>=0
                                                        SD(1)=0;
                                                        SD(2)=dataSD(x);
                                                    else
                                                        SD(1)=dataSD(x);
                                                        SD(2)=0;
                                                    end
                                                    errorbar(x,dataMeans(x),SD(1),SD(2),'k','LineStyle','none')
                                                end

                                                if indivLines==1
                                                    xl=1:numel(modalitiesRM{nRm});
                                                    plot(xl,data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nModInd1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nModInd2})(:,col4means{nRm1}(:,nModRm)),'k--'); hold on
                                                    scatter(xl,data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nModInd1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nModInd2})(:,col4means{nRm1}(:,nModRm)),'k+')
                                                end

                                                if abs(dataMeans(x))<1
                                                    text(x,0.5*(dataMeans(x)),sprintf('%.3f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                                elseif abs(dataMeans(x))<10
                                                    text(x,0.5*(dataMeans(x)),sprintf('%.2f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                                elseif abs(dataMeans(x))<100
                                                    text(x,0.5*(dataMeans(x)),sprintf('%.1f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                                else
                                                    text(x,0.5*(dataMeans(x)),sprintf('%.0f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                                end
                                            end

                                            xticks(1:numel(cond4effect{rmEffect(nRm1)}))
                                            xticklabels(cond4effect{rmEffect(nRm1)})
                                            if nModRm==numel(cond4effect{rmEffect(nRm2)})
                                                xlabel(condNames{rmEffect(nRm1)})
                                            end
                                            ylabel(units)
                                            box off
                                            ax{nModRm}=gca;
                                            ax{nModRm}.XGrid='off';
                                            ax{nModRm}.YGrid='on';
                                            yl(:,nModRm)=ylim;
                                            title(cond4effect{rmEffect(nRm2)}{nModRm})

                                        end

                                        nAovInt4=findComboStrings(aov.Effect,{condNames{indEffect(nInd1)}, condNames{indEffect(nInd2)}, condNames{rmEffect(nRm1)}, condNames{rmEffect(nRm2)}});
                                        pINT4=aov{nAovInt4,6};

                                        amp=[max(max(yl))-min(min(yl))];
                                        for nModRm=1:numel(cond4effect{rmEffect(nRm2)})

                                            dataMeans=nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nModInd1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nModInd2})(:,col4means{nRm1}(:,nModRm)));
                                            dataSD=nanstd(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nModInd1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nModInd2})(:,col4means{nRm1}(:,nModRm)));

                                            dataMeans4pv=dataMeans;
                                            if plotSD==1
                                                dataMeans4pv=dataMeans+sign(dataMeans).*dataSD;
                                            end
                                            if indivLines==1
                                                dataMeans4pv=max(abs(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nModInd1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nModInd2})(:,col4means{nRm1}(:,nModRm)))).*sign(dataMeans);
                                            end
                                            if plotSD==1 && indivLines==1
                                                dataMeans4pv=sign(dataMeans).*max([abs(dataMeans+sign(dataMeans).*dataSD); max(abs(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nModInd1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nModInd2})(:,col4means{nRm1}(:,nModRm)))).*sign(dataMeans)]);
                                            end

                                            isSignificant=0;
                                            set(f,'CurrentAxes',ax{nModRm});
                                            pValues=ones(1,size(postHoc.(condNamesVerif{rmEffect(nRm1)}),1));
                                            pSelected=0;
                                            if pINT4<pcritical(1)
                                                isSignificant=1;
                                                ph4=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]);
                                                phCut1=findcol(ph4{:,1}, allMod.(condNamesVerif{indEffect(nInd1)}){nModInd1});
                                                phCut2=findcol(ph4{:,2}, allMod.(condNamesVerif{indEffect(nInd2)}){nModInd2});
                                                phCut3=findcol(ph4{:,3}, allMod.(condNamesVerif{rmEffect(nRm2)}){nModRm});
                                                phCut4=intersect(phCut1, phCut2);
                                                phCut=intersect(phCut4, phCut3);
                                                pValues=ph4{phCut,8};
                                                pSelected=1;
                                            end
                                            if pINT4<pcritical(end) & pSelected==0
                                                ph4=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]);
                                                phCut1=findcol(ph4{:,1}, allMod.(condNamesVerif{indEffect(nInd1)}){nModInd1});
                                                phCut2=findcol(ph4{:,2}, allMod.(condNamesVerif{indEffect(nInd2)}){nModInd2});
                                                phCut3=findcol(ph4{:,3}, allMod.(condNamesVerif{rmEffect(nRm2)}){nModRm});
                                                phCut4=intersect(phCut1, phCut2);
                                                phCut=intersect(phCut4, phCut3);
                                                pValues=ph4{phCut,8};
                                                pSelected=1;
                                            end

                                            nSignificant=1;
                                            nColSignificant=ones(1,numel(dataMeans));
                                            if pSelected==1
                                                for i=1:numel(pValues)
                                                    pV=pValues(i);
                                                    if pV<=pcritical(end)
                                                        if abs(yl(2))>abs(yl(1))
                                                            if abs(dataMeans(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1)))>abs(dataMeans(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2)))
                                                                addPvalue(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1), dataMeans4pv(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1))+0.0225*numel(cond4effect{rmEffect(nRm2)})*nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1))*amp, pV, pcritical, colors{rmEffect(nRm1)}(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2),:))
                                                                nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1))=nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1))+1;
                                                            else
                                                                addPvalue(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2), dataMeans4pv(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2))+0.0225*numel(cond4effect{rmEffect(nRm2)})*nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2))*amp, pV, pcritical, colors{rmEffect(nRm1)}(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1),:))
                                                                nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2))=nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2))+1;
                                                            end
                                                        else
                                                            if abs(dataMeans(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1)))>abs(dataMeans(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2)))
                                                                addPvalue(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1), dataMeans4pv(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1))-0.0225*numel(cond4effect{rmEffect(nRm2)})*nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1))*amp, pV, pcritical, colors{rmEffect(nRm1)}(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2),:))
                                                                nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1))=nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1))+1;
                                                            else
                                                                addPvalue(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2), dataMeans4pv(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2))-0.0225*numel(cond4effect{rmEffect(nRm2)})*nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2))*amp, pV, pcritical, colors{rmEffect(nRm1)}(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,1),:))
                                                                nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2))=nColSignificant(order4ES.(condNamesVerif{rmEffect(nRm1)})(i,2))+1;
                                                            end
                                                        end
                                                        nSignificant=nSignificant+1;
                                                    end
                                                end
                                            end

                                            yl(:,nModRm)=ylim*1.05.^max(nColSignificant);

                                        end

                                        for nModRm=1:numel(ax)
                                            set(ax{nModRm},'ylim',[min(min(yl)) max(max(yl))])
                                        end

                                        print('-dtiff',['-r' num2str(imageResolution)],fullfile(saveDir, 'Text', condNames{rmEffect(nRm1)}, [condNames{indEffect(nInd1)} ' By ' condNames{indEffect(nInd2)} ' = ' allMod.(condNamesVerif{indEffect(nInd1)}){nModInd1} ' ' allMod.(condNamesVerif{indEffect(nInd2)}){nModInd2} ' By ' condNames{rmEffect(nRm2)}]))
                                        close
                                        clear ax yl dataMeans dataSD

                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    %% independant effects
    %% MAIN EFFECT
    if numel(indEffect)>0
        % Lines
        if statsLines
            for nInd=1:numel(indEffect)

                f=figure('units','centimeters','position',[0 0 6+4*numel(modalitiesInd{nInd}) 4+9/16*4*numel(modalitiesInd{nInd})],'visible','off');

                for nMod=1:numel(modalitiesInd{nInd})
                    dataMeans(nMod)=nanmean(nanmean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod})),2));
                    dataSD(nMod)=nanstd(nanmean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod})),2));
                end

                for x=1:numel(modalitiesInd{nInd})

                    h=bar(x,dataMeans(x)); hold on
                    h(1).FaceColor=colors{indEffect(nInd)}(x,:);

                    if plotSD==1
                        if dataMeans(x)>=0
                            SD(1)=0;
                            SD(2)=dataSD(x);
                        else
                            SD(1)=dataSD(x);
                            SD(2)=0;
                        end
                        errorbar(x,dataMeans(x),SD(1),SD(2),'k','LineStyle','none')
                    end

                    if indivLines==1
                        scatter(x,mean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){x})),2), 'kx' ,'handlevisibility','off')
                    end

                    if abs(dataMeans(x))<1
                        text(x,0.5*(dataMeans(x)),sprintf('%.3f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                    elseif abs(dataMeans(x))<10
                        text(x,0.5*(dataMeans(x)),sprintf('%.2f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                    elseif abs(dataMeans(x))<100
                        text(x,0.5*(dataMeans(x)),sprintf('%.1f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                    else
                        text(x,0.5*(dataMeans(x)),sprintf('%.0f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                    end

                end
                xticks(1:numel(modalitiesInd{nInd}))
                xticklabels(modalitiesInd{nInd})
                xlabel(condNames{indEffect(nInd)})
                ylabel(units)
                yl=ylim;
                box off
                ax=gca;
                ax.XGrid='off';
                ax.YGrid='on';

                nAov=findcolExact(aov.Effect, verifFieldName(condNames{indEffect(nInd)}));
                pMAIN=aov{nAov,6};

                amp=[max(yl)-min(yl)];
                isSignificant=0;
                pValues=ones(1,size(postHoc.(condNamesVerif{indEffect(nInd)}),1));
                pSelected=0;
                if pMAIN<pcritical(1) &  pSelected==0
                    isSignificant=1;
                    pValues=postHoc.(condNamesVerif{indEffect(nInd)}){:,5};
                    pSelected=1;
                end

                if pMAIN<pcritical(end) &  pSelected==0
                    pValues=postHoc.(condNamesVerif{indEffect(nInd)}){:,5};
                    pSelected=1;
                end

                nSignificant=1;
                if pSelected==1
                    nSignificant=1;
                    for i=1:numel(pValues)
                        pV=pValues(i);
                        if pV<=pcritical(end)
                            if abs(yl(2))>abs(yl(1))
                                if pV<pcritical(1) & isSignificant==1
                                    hline(yl(2)+0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNamesVerif{indEffect(nInd)})(i,:),'lineWidth',1.5);
                                else
                                    hline(yl(2)+0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNamesVerif{indEffect(nInd)})(i,:),'lineWidth',1.5);
                                end
                            else
                                if pV<pcritical(1) & isSignificant==1
                                    hline(yl(1)-0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNamesVerif{indEffect(nInd)})(i,:),'lineWidth',1.5);
                                else
                                    hline(yl(1)-0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNamesVerif{indEffect(nInd)})(i,:),'lineWidth',1.5);
                                end
                            end
                            nSignificant=nSignificant+1;
                        end
                    end
                end


                yl=yl*1.05.^max(nSignificant);
                set(ax,'ylim',[min(yl) max(yl)])

                print('-dtiff',['-r' num2str(imageResolution)], fullfile(saveDir, 'Lines', condNames{indEffect(nInd)}, 'All RM'))
                close


                if indivLines==2
                    f=figure('units','centimeters','position',[0 0 6+4*numel(modalitiesInd{nInd}) 4+9/16*4*numel(modalitiesInd{nInd})],'visible','off');
                    hline(0, 'xLimits',[0 (numel(modalitiesRM{nRm})+1)],'linetype','-k', 'linewidth',0.5); hold on
                    for nMod=1:numel(allMod.(condNamesVerif{indEffect(nInd)}))
                        for n=1:size(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod})),1)
                            if n==1
                                scatter(nMod,mean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(n,:)), '+', 'MarkerFaceColor', colors{indEffect(nInd)}(nMod,:), 'MarkerEdgeColor', colors{indEffect(nInd)}(nMod,:),'handlevisibility','on')
                            else
                                scatter(nMod,mean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(n,:)), '+', 'MarkerFaceColor', colors{indEffect(nInd)}(nMod,:), 'MarkerEdgeColor', colors{indEffect(nInd)}(nMod,:),'handlevisibility','off')
                            end
                        end
                    end
                    xlim(xlp)
                    xticks(1:numel(allMod.(condNamesVerif{indEffect(nInd)})))
                    xticklabels(allMod.(condNamesVerif{indEffect(nInd)}))
                    xlabel(condNames{indEffect(nInd)})
                    ylabel(units)
                    yl=ylim;
                    ylim([-max(abs(yl)) max(abs(yl))])
                    box off
                    ax=gca;
                    ax.XGrid='off';
                    ax.YGrid='on';
                    legend(allModalities{indEffect(nInd)}, "box", "off", "Location", "northeast");

                    print('-dtiff',['-r' num2str(imageResolution)],fullfile(saveDir, 'Lines', condNames{indEffect(nInd)}, 'All RM and indiv'))
                    close
                end

                clear ax yl dataMeans dataSD

            end
        end

        % Text
        for nInd=1:numel(indEffect)

            f=figure('units','centimeters','position',[0 0 6+4*numel(modalitiesInd{nInd}) 4+9/16*4*numel(modalitiesInd{nInd})],'visible','off');

            for nMod=1:numel(modalitiesInd{nInd})
                dataMeans(nMod)=nanmean(nanmean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod})),2));
                dataSD(nMod)=nanstd(nanmean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod})),2));
            end

            for x=1:numel(modalitiesInd{nInd})

                h=bar(x,dataMeans(x)); hold on
                h(1).FaceColor=colors{indEffect(nInd)}(x,:);

                if plotSD==1
                    if dataMeans(x)>=0
                        SD(1)=0;
                        SD(2)=dataSD(x);
                    else
                        SD(1)=dataSD(x);
                        SD(2)=0;
                    end
                    errorbar(x,dataMeans(x),SD(1),SD(2),'k','LineStyle','none')
                end

                if indivLines==1
                    scatter(x,mean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){x})),2), 'kx' ,'handlevisibility','off')
                end

                if abs(dataMeans(x))<1
                    text(x,0.5*(dataMeans(x)),sprintf('%.3f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                elseif abs(dataMeans(x))<10
                    text(x,0.5*(dataMeans(x)),sprintf('%.2f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                elseif abs(dataMeans(x))<100
                    text(x,0.5*(dataMeans(x)),sprintf('%.1f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                else
                    text(x,0.5*(dataMeans(x)),sprintf('%.0f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                end

            end
            xticks(1:numel(modalitiesInd{nInd}))
            xticklabels(modalitiesInd{nInd})
            xlabel(condNames{indEffect(nInd)})
            ylabel(units)
            yl=ylim;
            box off
            ax=gca;
            ax.XGrid='off';
            ax.YGrid='on';

            if plotSD==1
                dataMeans=dataMeans+sign(dataMeans).*dataSD;
            end

            nAov=findcolExact(aov.Effect, verifFieldName(condNames{indEffect(nInd)}));
            pMAIN=aov{nAov,6};
            if iscell(pMAIN)
                pMAIN=pMAIN{1};
            end
            amp=[max(yl)-min(yl)];
            isSignificant=0;
            pValues=ones(1,size(postHoc.(condNamesVerif{indEffect(nInd)}),1));
            pSelected=0;
            if pMAIN<pcritical(1) &  pSelected==0
                isSignificant=1;
                pValues=postHoc.(condNamesVerif{indEffect(nInd)}){:,5};
                pSelected=1;
            end

            if pMAIN<pcritical(end) &  pSelected==0
                pValues=postHoc.(condNamesVerif{indEffect(nInd)}){:,5};
                pSelected=1;
            end

            nColSignificant=ones(1,numel(dataMeans));
            nSignificant=1;
            if pSelected==1
                for i=1:numel(pValues)
                    pV=pValues(i);
                    if pV<=pcritical(end)
                        if abs(yl(2))>abs(yl(1))
                            if abs(dataMeans(order4ES.(condNamesVerif{indEffect(nInd)})(i,1)))>abs(dataMeans(order4ES.(condNamesVerif{indEffect(nInd)})(i,2)))
                                addPvalue(order4ES.(condNamesVerif{indEffect(nInd)})(i,1), yl(2)+0.055*nColSignificant(order4ES.(condNamesVerif{indEffect(nInd)})(i,1))*amp, pV, pcritical, colors{indEffect(nInd)}(order4ES.(condNamesVerif{indEffect(nInd)})(i,2),:))
                                nColSignificant(order4ES.(condNamesVerif{indEffect(nInd)})(i,1))=nColSignificant(order4ES.(condNamesVerif{indEffect(nInd)})(i,1))+1;
                            else
                                addPvalue(order4ES.(condNamesVerif{indEffect(nInd)})(i,2), yl(2)+0.055*nColSignificant(order4ES.(condNamesVerif{indEffect(nInd)})(i,2))*amp, pV, pcritical, colors{indEffect(nInd)}(order4ES.(condNamesVerif{indEffect(nInd)})(i,1),:))
                                nColSignificant(order4ES.(condNamesVerif{indEffect(nInd)})(i,2))=nColSignificant(order4ES.(condNamesVerif{indEffect(nInd)})(i,2))+1;
                            end
                        else
                            if abs(dataMeans(order4ES.(condNamesVerif{indEffect(nInd)})(i,1)))>abs(dataMeans(order4ES.(condNamesVerif{indEffect(nInd)})(i,2)))
                                addPvalue(order4ES.(condNamesVerif{indEffect(nInd)})(i,1), yl(1)-0.055*nColSignificant(order4ES.(condNamesVerif{indEffect(nInd)})(i,1))*amp, pV, pcritical, colors{indEffect(nInd)}(order4ES.(condNamesVerif{indEffect(nInd)})(i,2),:))
                                nColSignificant(order4ES.(condNamesVerif{indEffect(nInd)})(i,1))=nColSignificant(order4ES.(condNamesVerif{indEffect(nInd)})(i,1))+1;
                            else
                                addPvalue(order4ES.(condNamesVerif{indEffect(nInd)})(i,2), yl(1)-0.055*nColSignificant(order4ES.(condNamesVerif{indEffect(nInd)})(i,2))*amp, pV, pcritical, colors{indEffect(nInd)}(order4ES.(condNamesVerif{indEffect(nInd)})(i,1),:))
                                nColSignificant(order4ES.(condNamesVerif{indEffect(nInd)})(i,2))=nColSignificant(order4ES.(condNamesVerif{indEffect(nInd)})(i,2))+1;
                            end
                        end
                        nSignificant=nSignificant+1;
                    end
                end
            end

            yl=yl*1.05.^max(nColSignificant);
            set(ax,'ylim',[min(yl) max(yl)])
            xlp=xlim;

            print('-dtiff',['-r' num2str(imageResolution)], fullfile(saveDir, 'text', condNames{indEffect(nInd)}, 'All RM'))
            close

            if indivLines==2
                f=figure('units','centimeters','position',[0 0 6+4*numel(modalitiesInd{nInd}) 4+9/16*4*numel(modalitiesInd{nInd})],'visible','off');
                hline(0, 'xLimits',[0 (numel(modalitiesRM{nRm})+1)],'linetype','-k', 'linewidth',0.5); hold on
                for nMod=1:numel(allMod.(condNamesVerif{indEffect(nInd)}))
                    for n=1:size(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod})),1)
                        if n==1
                            scatter(nMod,mean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(n,:)), '+', 'MarkerFaceColor', colors{indEffect(nInd)}(nMod,:), 'MarkerEdgeColor', colors{indEffect(nInd)}(nMod,:),'handlevisibility','on')
                        else
                            scatter(nMod,mean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){nMod}))(n,:)), '+', 'MarkerFaceColor', colors{indEffect(nInd)}(nMod,:), 'MarkerEdgeColor', colors{indEffect(nInd)}(nMod,:),'handlevisibility','off')
                        end
                    end
                end
                xlim(xlp)
                xticks(1:numel(allMod.(condNamesVerif{indEffect(nInd)})))
                xticklabels(allMod.(condNamesVerif{indEffect(nInd)}))
                xlabel(condNames{indEffect(nInd)})
                ylabel(units)
                yl=ylim;
                ylim([-max(abs(yl)) max(abs(yl))])
                box off
                ax=gca;
                ax.XGrid='off';
                ax.YGrid='on';
                legend(allModalities{indEffect(nInd)}, "box", "off", "Location", "northeast");

                print('-dtiff',['-r' num2str(imageResolution)],fullfile(saveDir, 'Text', condNames{indEffect(nInd)}, 'All RM and indiv'))
                close
            end
        end
        clear ax yl dataMeans dataSD

    end

    % 1 IND 1RM
    if numel(indEffect)>0 & numel(rmEffect)>0
        % Lines
        if statsLines
            for nRm=1:numel(effectRM)
                for nInd=1:numel(indEffect)

                    f=figure('units','centimeters','position',[0 0 6+4*numel(modalitiesInd{nInd}) 4+9/16*4*numel(modalitiesInd{nInd})*numel(modalitiesRM{nRm})],'visible','off');

                    for nMod=1:numel(allMod.(condNamesVerif{rmEffect(nRm)}))

                        subplot(numel(allMod.(condNamesVerif{rmEffect(nRm)})),1,nMod)

                        for x=1:numel(modalitiesInd{nInd})
                            dataMeans(x)=nanmean(nanmean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){x}))(:,col4means{nRm}(nMod,:)),2));
                            dataSD(x)=nanstd(nanmean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){x}))(:,col4means{nRm}(nMod,:)),2));
                        end

                        for x=1:numel(modalitiesInd{nInd})

                            h=bar(x,dataMeans(x)); hold on
                            h(1).FaceColor=colors{indEffect(nInd)}(x,:);

                            if plotSD==1
                                if dataMeans(x)>=0
                                    SD(1)=0;
                                    SD(2)=dataSD(x);
                                else
                                    SD(1)=dataSD(x);
                                    SD(2)=0;
                                end
                                errorbar(x,dataMeans(x),SD(1),SD(2),'k','LineStyle','none')
                            end

                            if abs(dataMeans(x))<1
                                text(x,0.5*(dataMeans(x)),sprintf('%.3f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                            elseif abs(dataMeans(x))<10
                                text(x,0.5*(dataMeans(x)),sprintf('%.2f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                            elseif abs(dataMeans(x))<100
                                text(x,0.5*(dataMeans(x)),sprintf('%.1f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                            else
                                text(x,0.5*(dataMeans(x)),sprintf('%.0f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                            end

                        end
                        xticks(1:numel(modalitiesInd{nInd}))
                        xticklabels(modalitiesInd{nInd})
                        if nMod==numel(allMod.(condNamesVerif{rmEffect(nRm)}))
                            xlabel(condNames{indEffect(nInd)})
                        end
                        ylabel(units)
                        yl(:,nMod)=ylim;
                        box off
                        ax{nMod}=gca;
                        ax{nMod}.XGrid='off';
                        ax{nMod}.YGrid='on';
                        title(allMod.(condNamesVerif{rmEffect(nRm)}){nMod})
                    end

                    nAovInt=findcolExact(aov.Effect,[condNames{indEffect(nInd)} ':' condNames{rmEffect(nRm)}]);
                    nAov=findcolExact(aov.Effect, verifFieldName(condNames{indEffect(nInd)}));
                    pINT=aov{nAovInt,6};
                    pMAIN=aov{nAov,6};

                    amp=[max(max(yl))-min(min(yl))];

                    for nMod=1:numel(allMod.(condNamesVerif{rmEffect(nRm)}))
                        isSignificant=0;
                        set(f,'CurrentAxes',ax{nMod});

                        for x=1:numel(modalitiesInd{nInd})
                            dataMeans(x)=nanmean(nanmean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){x}))(:,col4means{nRm}(nMod,:)),2));
                            dataSD(x)=nanstd(nanmean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){x}))(:,col4means{nRm}(nMod,:)),2));
                        end

                        pValues=ones(1,size(postHoc.(condNamesVerif{indEffect(nInd)}),1));
                        pSelected=0;
                        if pINT<pcritical(1)
                            isSignificant=1;
                            rows=findcolExact(postHoc.([condNames{indEffect(nInd)} 'By' condNames{rmEffect(nRm)}]){:,1},allMod.(condNamesVerif{rmEffect(nRm)}){nMod});
                            pValues=postHoc.([condNames{indEffect(nInd)} 'By' condNames{rmEffect(nRm)}]){rows,6};
                            pSelected=1;
                        end
                        %                     if pMAIN<pcritical(1) &  pSelected==0
                        %                         isSignificant=1;
                        %                         pValues=postHoc.(condNamesVerif{indEffect(nInd)}){:,5};
                        %                         pSelected=1;
                        %                     end
                        if pINT<pcritical(end) &  pSelected==0
                            rows=findcolExact(postHoc.([condNames{indEffect(nInd)} 'By' condNames{rmEffect(nRm)}]){:,1},allMod.(condNamesVerif{rmEffect(nRm)}){nMod});
                            pValues=postHoc.([condNames{indEffect(nInd)} 'By' condNames{rmEffect(nRm)}]){rows,6};
                            pSelected=1;
                        end
                        %                     if pMAIN<pcritical(end) &  pSelected==0
                        %                         pValues=postHoc.(condNamesVerif{indEffect(nInd)}){:,5};
                        %                         pSelected=1;
                        %                     end

                        if plotSD==1
                            dataMeans=dataMeans+sign(dataMeans).*dataSD;
                        end

                        nSignificant=1;
                        if pSelected==1
                            nSignificant=1;
                            for i=1:numel(pValues)
                                pV=pValues(i);
                                if pV<=pcritical(end)
                                    if abs(yl(2))>abs(yl(1))
                                        if pV<pcritical(1) & isSignificant==1
                                            hline(yl(2)+0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNamesVerif{indEffect(nInd)})(i,:),'lineWidth',1.5);
                                        else
                                            hline(yl(2)+0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNamesVerif{indEffect(nInd)})(i,:),'lineWidth',1.5);
                                        end
                                    else
                                        if pV<pcritical(1) & isSignificant==1
                                            hline(yl(1)-0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNamesVerif{indEffect(nInd)})(i,:),'lineWidth',1.5);
                                        else
                                            hline(yl(1)-0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNamesVerif{indEffect(nInd)})(i,:),'lineWidth',1.5);
                                        end
                                    end
                                    nSignificant=nSignificant+1;
                                end
                            end
                        end
                    end

                    yl=yl*1.05.^max(nSignificant);
                    for nMod=1:numel(allMod.(condNamesVerif{rmEffect(nRm)}))
                        set(ax{nMod},'ylim',[min(min(yl)) max(max(yl))])
                    end

                    print('-dtiff',['-r' num2str(imageResolution)], fullfile(saveDir, 'Lines', condNames{indEffect(nInd)}, ['By ' condNames{rmEffect(nRm)}]))
                    close
                    clear ax yl dataMeans dataSD

                end
            end
        end

        % Text
        for nRm=1:numel(effectRM)
            for nInd=1:numel(indEffect)

                f=figure('units','centimeters','position',[0 0 6+4*numel(modalitiesInd{nInd}) 4+9/16*4*numel(modalitiesInd{nInd})*numel(modalitiesRM{nRm})],'visible','off');

                for nMod=1:numel(allMod.(condNamesVerif{rmEffect(nRm)}))

                    subplot(numel(allMod.(condNamesVerif{rmEffect(nRm)})),1,nMod)

                    for x=1:numel(modalitiesInd{nInd})
                        dataMeans(x)=nanmean(nanmean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){x}))(:,col4means{nRm}(nMod,:)),2));
                        dataSD(x)=nanstd(nanmean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){x}))(:,col4means{nRm}(nMod,:)),2));
                    end

                    for x=1:numel(modalitiesInd{nInd})

                        h=bar(x,dataMeans(x)); hold on
                        h(1).FaceColor=colors{indEffect(nInd)}(x,:);

                        if plotSD==1
                            if dataMeans(x)>=0
                                SD(1)=0;
                                SD(2)=dataSD(x);
                            else
                                SD(1)=dataSD(x);
                                SD(2)=0;
                            end
                            errorbar(x,dataMeans(x),SD(1),SD(2),'k','LineStyle','none')
                        end

                        if abs(dataMeans(x))<1
                            text(x,0.5*(dataMeans(x)),sprintf('%.3f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                        elseif abs(dataMeans(x))<10
                            text(x,0.5*(dataMeans(x)),sprintf('%.2f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                        elseif abs(dataMeans(x))<100
                            text(x,0.5*(dataMeans(x)),sprintf('%.1f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                        else
                            text(x,0.5*(dataMeans(x)),sprintf('%.0f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                        end

                    end
                    xticks(1:numel(modalitiesInd{nInd}))
                    xticklabels(modalitiesInd{nInd})
                    if nMod==numel(allMod.(condNamesVerif{rmEffect(nRm)}))
                        xlabel(condNames{indEffect(nInd)})
                    end
                    ylabel(units)
                    yl(:,nMod)=ylim;
                    box off
                    ax{nMod}=gca;
                    ax{nMod}.XGrid='off';
                    ax{nMod}.YGrid='on';
                    title(allMod.(condNamesVerif{rmEffect(nRm)}){nMod})
                end

                nAovInt=findcolExact(aov.Effect,[condNames{indEffect(nInd)} ':' condNames{rmEffect(nRm)}]);
                nAov=findcolExact(aov.Effect, verifFieldName(condNames{indEffect(nInd)}));
                pINT=aov{nAovInt,6};
                pMAIN=aov{nAov,6};

                amp=[max(max(yl))-min(min(yl))];

                for nMod=1:numel(allMod.(condNamesVerif{rmEffect(nRm)}))

                    isSignificant=0;
                    set(f,'CurrentAxes',ax{nMod});

                    for x=1:numel(modalitiesInd{nInd})
                        dataMeans(x)=nanmean(nanmean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){x}))(:,col4means{nRm}(nMod,:)),2));
                        dataSD(x)=nanstd(nanmean(data4plot.(condNamesVerif{indEffect(nInd)}).(verifFieldName(allMod.(condNamesVerif{indEffect(nInd)}){x}))(:,col4means{nRm}(nMod,:)),2));
                    end

                    pValues=ones(1,size(postHoc.(condNamesVerif{indEffect(nInd)}),1));
                    pSelected=0;
                    if pINT<pcritical(1)
                        isSignificant=1;
                        rows=findcolExact(postHoc.([condNames{indEffect(nInd)} 'By' condNames{rmEffect(nRm)}]){:,1},allMod.(condNamesVerif{rmEffect(nRm)}){nMod});
                        pValues=postHoc.([condNames{indEffect(nInd)} 'By' condNames{rmEffect(nRm)}]){rows,6};
                        pSelected=1;
                    end
                    %                 if pMAIN<pcritical(1) &  pSelected==0
                    %                     isSignificant=1;
                    %                     pValues=postHoc.(condNamesVerif{indEffect(nInd)}){:,5};
                    %                     pSelected=1;
                    %                 end
                    if pINT<pcritical(end) &  pSelected==0
                        rows=findcolExact(postHoc.([condNames{indEffect(nInd)} 'By' condNames{rmEffect(nRm)}]){:,1},allMod.(condNamesVerif{rmEffect(nRm)}){nMod});
                        pValues=postHoc.([condNames{indEffect(nInd)} 'By' condNames{rmEffect(nRm)}]){rows,6};
                        pSelected=1;
                    end
                    %                 if pMAIN<pcritical(end) &  pSelected==0
                    %                     pValues=postHoc.(condNamesVerif{indEffect(nInd)}){:,5};
                    %                     pSelected=1;
                    %                 end

                    dataMeans4pv=dataMeans;
                    if plotSD==1
                        dataMeans4pv=dataMeans+sign(dataMeans).*dataSD;
                    end

                    nColSignificant=ones(1,numel(dataMeans));
                    nSignificant=1;
                    if pSelected==1
                        for i=1:numel(pValues)
                            pV=pValues(i);
                            if pV<=pcritical(end)
                                if yl(2,nMod)>0
                                    if (dataMeans(order4ES.(condNamesVerif{indEffect(nInd)})(i,1)))>(dataMeans(order4ES.(condNamesVerif{indEffect(nInd)})(i,2)))
                                        addPvalue(order4ES.(condNamesVerif{indEffect(nInd)})(i,1), dataMeans4pv(order4ES.(condNamesVerif{indEffect(nInd)})(i,1))+0.0225*numel(allMod.(condNamesVerif{rmEffect(nRm)}))*nColSignificant(order4ES.(condNamesVerif{indEffect(nInd)})(i,1))*amp, pV, pcritical, colors{indEffect(nInd)}(order4ES.(condNamesVerif{indEffect(nInd)})(i,2),:))
                                        nColSignificant(order4ES.(condNamesVerif{indEffect(nInd)})(i,1))=nColSignificant(order4ES.(condNamesVerif{indEffect(nInd)})(i,1))+1;
                                    else
                                        addPvalue(order4ES.(condNamesVerif{indEffect(nInd)})(i,2), dataMeans4pv(order4ES.(condNamesVerif{indEffect(nInd)})(i,2))+0.0225*numel(allMod.(condNamesVerif{rmEffect(nRm)}))*nColSignificant(order4ES.(condNamesVerif{indEffect(nInd)})(i,2))*amp, pV, pcritical, colors{indEffect(nInd)}(order4ES.(condNamesVerif{indEffect(nInd)})(i,1),:))
                                        nColSignificant(order4ES.(condNamesVerif{indEffect(nInd)})(i,2))=nColSignificant(order4ES.(condNamesVerif{indEffect(nInd)})(i,2))+1;
                                    end
                                else
                                    if abs(dataMeans(order4ES.(condNamesVerif{indEffect(nInd)})(i,1)))>abs(dataMeans(order4ES.(condNamesVerif{indEffect(nInd)})(i,2)))
                                        addPvalue(order4ES.(condNamesVerif{indEffect(nInd)})(i,1), dataMeans4pv(order4ES.(condNamesVerif{indEffect(nInd)})(i,1))-0.0225*numel(allMod.(condNamesVerif{rmEffect(nRm)}))*nColSignificant(order4ES.(condNamesVerif{indEffect(nInd)})(i,1))*amp, pV, pcritical, colors{indEffect(nInd)}(order4ES.(condNamesVerif{indEffect(nInd)})(i,2),:))
                                        nColSignificant(order4ES.(condNamesVerif{indEffect(nInd)})(i,1))=nColSignificant(order4ES.(condNamesVerif{indEffect(nInd)})(i,1))+1;
                                    else
                                        addPvalue(order4ES.(condNamesVerif{indEffect(nInd)})(i,2), dataMeans4pv(order4ES.(condNamesVerif{indEffect(nInd)})(i,2))-0.0225*numel(allMod.(condNamesVerif{rmEffect(nRm)}))*nColSignificant(order4ES.(condNamesVerif{indEffect(nInd)})(i,2))*amp, pV, pcritical, colors{indEffect(nInd)}(order4ES.(condNamesVerif{indEffect(nInd)})(i,1),:))
                                        nColSignificant(order4ES.(condNamesVerif{indEffect(nInd)})(i,2))=nColSignificant(order4ES.(condNamesVerif{indEffect(nInd)})(i,2))+1;
                                    end
                                end
                                nSignificant=nSignificant+1;
                            end
                        end
                    end
                end

                yl=yl*1.05.^max(nColSignificant);
                for nMod=1:numel(allMod.(condNamesVerif{rmEffect(nRm)}))
                    set(ax{nMod},'ylim',[min(min(yl)) max(max(yl))])
                end

                print('-dtiff',['-r' num2str(imageResolution)], fullfile(saveDir, 'Text', condNames{indEffect(nInd)}, ['By ' condNames{rmEffect(nRm)}]))
                close
                clear ax yl dataMeans dataSD

            end
        end
    end

    %% Interaction of 2 IND effect
    if numel(indEffect)>1
        % Lines
        if statsLines
            for nInd1=1:numel(indEffect)
                for nInd2=1:numel(indEffect)
                    if nInd1~=nInd2

                        f=figure('units','centimeters','position',[0 0 6+4*numel(modalitiesInd{nInd1}) 4+9/16*4*numel(modalitiesInd{nInd1})*numel(modalitiesInd{nInd2})],'visible','off');

                        for nMod2=1:numel(modalitiesInd{nInd2})

                            subplot(numel(modalitiesInd{nInd2}),1,nMod2)

                            for nMod1=1:numel(modalitiesInd{nInd1})
                                dataMeans(nMod1)=nanmean(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nMod1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2}),2));
                                dataSD(nMod1)=nanstd(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nMod1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2}),2));
                            end

                            for x=1:numel(modalitiesInd{nInd1})

                                h=bar(x,dataMeans(x)); hold on
                                h(1).FaceColor=colors{indEffect(nInd1)}(x,:);

                                if plotSD==1
                                    if dataMeans(x)>=0
                                        SD(1)=0;
                                        SD(2)=dataSD(x);
                                    else
                                        SD(1)=dataSD(x);
                                        SD(2)=0;
                                    end
                                    errorbar(x,dataMeans(x),SD(1),SD(2),'k','LineStyle','none')
                                end

                                if abs(dataMeans(x))<1
                                    text(x,0.5*(dataMeans(x)),sprintf('%.3f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                elseif abs(dataMeans(x))<10
                                    text(x,0.5*(dataMeans(x)),sprintf('%.2f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                elseif abs(dataMeans(x))<100
                                    text(x,0.5*(dataMeans(x)),sprintf('%.1f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                else
                                    text(x,0.5*(dataMeans(x)),sprintf('%.0f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                end

                            end
                            xticks(1:numel(modalitiesInd{nInd1}))
                            xticklabels(modalitiesInd{nInd1})
                            if nMod2==numel(modalitiesInd{nInd2})
                                xlabel(condNames{indEffect(nInd1)})
                            end
                            ylabel(units)
                            yl(:,nMod2)=ylim;
                            box off
                            ax{nMod2}=gca;
                            ax{nMod2}.XGrid='off';
                            ax{nMod2}.YGrid='on';
                            title(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2})

                        end

                        nAovInt=findcolExact(aov.Effect,[condNames{indEffect(nInd1)} ':' condNames{indEffect(nInd2)}]);
                        if isempty(nAovInt)
                            nAovInt=findcolExact(aov.Effect,[condNames{indEffect(nInd2)} ':' condNames{indEffect(nInd1)}]);
                        end
                        nAov=findcolExact(aov.Effect, verifFieldName(condNames{indEffect(nInd1)}));
                        pINT=aov{nAovInt,6};
                        pMAIN=aov{nAov,6};

                        amp=[max(max(yl))-min(min(yl))];

                        for nMod2=1:numel(modalitiesInd{nInd2})

                            for nMod1=1:numel(modalitiesInd{nInd1})
                                dataMeans(nMod1)=nanmean(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nMod1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2}),2));
                                dataSD(nMod1)=nanstd(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nMod1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2}),2));
                            end

                            isSignificant=0;
                            set(f,'CurrentAxes',ax{nMod2});
                            pValues=ones(1,size(postHoc.(condNamesVerif{indEffect(nInd1)}),1));
                            pSelected=0;
                            if pINT<pcritical(1)
                                isSignificant=1;
                                phCut=findcolExact(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,1}, allMod.(condNamesVerif{indEffect(nInd2)}){nMod2});
                                pValues=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){phCut,6};
                                pSelected=1;
                            end
                            %                         if pMAIN<pcritical(1) &  pSelected==0
                            %                             isSignificant=1;
                            %                             pValues=postHoc.(condNamesVerif{indEffect(nInd1)}){:,5};
                            %                             pSelected=1;
                            %                         end
                            if pINT<pcritical(end) &  pSelected==0
                                phCut=findcolExact(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,1}, allMod.(condNamesVerif{indEffect(nInd2)}){nMod2});
                                pValues=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){phCut,6};
                                pSelected=1;
                            end
                            %                         if pMAIN<pcritical(end) &  pSelected==0
                            %                             pValues=postHoc.(condNamesVerif{indEffect(nInd1)}){:,5};
                            %                             pSelected=1;
                            %                         end

                            if pSelected==1
                                nSignificant=1;
                                for i=1:numel(pValues)
                                    pV=pValues(i);
                                    if pV<=pcritical(end)
                                        if yl(2,nMod2)>0
                                            if pV<pcritical(1) & isSignificant==1
                                                hline(yl(2,nMod2)+0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNamesVerif{indEffect(nInd1)})(i,:),'lineWidth',1.5);
                                            else
                                                hline(yl(2,nMod2)+0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNamesVerif{indEffect(nInd1)})(i,:),'lineWidth',1.5);
                                            end
                                        else
                                            if pV<pcritical(1) & isSignificant==1
                                                hline(yl(1,nMod2)-0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNamesVerif{indEffect(nInd1)})(i,:),'lineWidth',1.5);
                                            else
                                                hline(yl(1,nMod2)-0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNamesVerif{indEffect(nInd1)})(i,:),'lineWidth',1.5);
                                            end
                                        end
                                        nSignificant=nSignificant+1;
                                    end
                                end
                            end

                            yl(:,nMod2)=ylim;

                        end

                        for nMod2=1:numel(ax)
                            set(ax{nMod2},'ylim',1.05*[min(min(yl)) max(max(yl))])
                        end

                        print('-dtiff',['-r' num2str(imageResolution)],fullfile(saveDir, 'Lines', condNames{indEffect(nInd1)}, ['All RM by ' condNames{indEffect(nInd2)}]))
                        close
                        clear ax yl dataMeans dataSD
                    end
                end
            end
        end

        % Text
        for nInd1=1:numel(indEffect)
            for nInd2=1:numel(indEffect)
                if nInd1~=nInd2

                    f=figure('units','centimeters','position',[0 0 6+4*numel(modalitiesInd{nInd1}) 4+9/16*4*numel(modalitiesInd{nInd1})*numel(modalitiesInd{nInd2})],'visible','off');

                    for nMod2=1:numel(modalitiesInd{nInd2})

                        subplot(numel(modalitiesInd{nInd2}),1,nMod2)

                        for nMod1=1:numel(modalitiesInd{nInd1})
                            dataMeans(nMod1)=nanmean(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nMod1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2}),2));
                            dataSD(nMod1)=nanstd(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nMod1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2}),2));
                        end

                        for x=1:numel(modalitiesInd{nInd1})

                            h=bar(x,dataMeans(x)); hold on
                            h(1).FaceColor=colors{indEffect(nInd1)}(x,:);

                            if plotSD==1
                                if dataMeans(x)>=0
                                    SD(1)=0;
                                    SD(2)=dataSD(x);
                                else
                                    SD(1)=dataSD(x);
                                    SD(2)=0;
                                end
                                errorbar(x,dataMeans(x),SD(1),SD(2),'k','LineStyle','none')
                            end

                            if abs(dataMeans(x))<1
                                text(x,0.5*(dataMeans(x)),sprintf('%.3f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                            elseif abs(dataMeans(x))<10
                                text(x,0.5*(dataMeans(x)),sprintf('%.2f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                            elseif abs(dataMeans(x))<100
                                text(x,0.5*(dataMeans(x)),sprintf('%.1f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                            else
                                text(x,0.5*(dataMeans(x)),sprintf('%.0f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                            end

                        end
                        xticks(1:numel(modalitiesInd{nInd1}))
                        xticklabels(modalitiesInd{nInd1})
                        if nMod2==numel(modalitiesInd{nInd2})
                            xlabel(condNames{indEffect(nInd1)})
                        end
                        ylabel(units)
                        yl(:,nMod2)=ylim;
                        box off
                        ax{nMod2}=gca;
                        ax{nMod2}.XGrid='off';
                        ax{nMod2}.YGrid='on';
                        title(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2})

                    end

                    nAovInt=findcolExact(aov.Effect,[condNames{indEffect(nInd1)} ':' condNames{indEffect(nInd2)}]);
                    if isempty(nAovInt)
                        nAovInt=findcolExact(aov.Effect,[condNames{indEffect(nInd2)} ':' condNames{indEffect(nInd1)}]);
                    end
                    nAov=findcolExact(aov.Effect, verifFieldName(condNames{indEffect(nInd1)}));
                    pINT=aov{nAovInt,6};
                    pMAIN=aov{nAov,6};

                    amp=[max(max(yl))-min(min(yl))];

                    for nMod2=1:numel(modalitiesInd{nInd2})

                        for nMod1=1:numel(modalitiesInd{nInd1})
                            dataMeans(nMod1)=nanmean(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nMod1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2}),2));
                            dataSD(nMod1)=nanstd(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nMod1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2}),2));
                        end

                        dataMeans4pv=dataMeans;
                        if plotSD==1
                            dataMeans4pv=dataMeans+sign(dataMeans).*dataSD;
                        end

                        isSignificant=0;
                        set(f,'CurrentAxes',ax{nMod2});
                        pValues=ones(1,size(postHoc.(condNamesVerif{indEffect(nInd1)}),1));
                        pSelected=0;
                        if pINT<pcritical(1)
                            isSignificant=1;
                            phCut=findcolExact(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,1}, allMod.(condNamesVerif{indEffect(nInd2)}){nMod2});
                            pValues=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){phCut,6};
                            pSelected=1;
                        end
                        %                     if pMAIN<pcritical(1) &  pSelected==0
                        %                         isSignificant=1;
                        %                         pValues=postHoc.(condNamesVerif{indEffect(nInd1)}){:,5};
                        %                         pSelected=1;
                        %                     end
                        if pINT<pcritical(end) &  pSelected==0
                            phCut=findcolExact(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,1}, allMod.(condNamesVerif{indEffect(nInd2)}){nMod2});
                            pValues=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){phCut,6};
                            pSelected=1;
                        end
                        %                     if pMAIN<pcritical(end) &  pSelected==0
                        %                         pValues=postHoc.(condNamesVerif{indEffect(nInd1)}){:,5};
                        %                         pSelected=1;
                        %                     end

                        nColSignificant=ones(1,numel(dataMeans));
                        nSignificant=1;
                        if pSelected==1
                            for i=1:numel(pValues)
                                pV=pValues(i);
                                if pV<=pcritical(end)
                                    if abs(yl(2))>abs(yl(1))
                                        if abs(dataMeans(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1)))>abs(dataMeans(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2)))
                                            addPvalue(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1), dataMeans4pv(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1))+0.0225*numel(modalitiesInd{nInd2})*nColSignificant(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1))*amp, pV, pcritical, colors{indEffect(nInd1)}(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2),:))
                                            nColSignificant(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1))=nColSignificant(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1))+1;
                                        else
                                            addPvalue(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2), dataMeans4pv(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2))+0.0225*numel(modalitiesInd{nInd2})*nColSignificant(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2))*amp, pV, pcritical, colors{indEffect(nInd1)}(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1),:))
                                            nColSignificant(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2))=nColSignificant(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2))+1;
                                        end
                                    else
                                        if abs(dataMeans(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1)))>abs(dataMeans(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2)))
                                            addPvalue(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1), dataMeans4pv(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1))-0.0225*numel(modalitiesInd{nInd2})*nColSignificant(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1))*amp, pV, pcritical, colors{indEffect(nInd1)}(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2),:))
                                            nColSignificant(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1))=nColSignificant(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1))+1;
                                        else
                                            addPvalue(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2), dataMeans4pv(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2))-0.0225*numel(modalitiesInd{nInd2})*nColSignificant(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2))*amp, pV, pcritical, colors{indEffect(nInd1)}(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1),:))
                                            nColSignificant(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2))=nColSignificant(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2))+1;
                                        end
                                    end
                                    nSignificant=nSignificant+1;
                                end
                            end
                        end
                    end

                    yl=yl*1.05.^max(nColSignificant);
                    for nMod2=1:numel(ax)
                        set(ax{nMod2},'ylim',1.05*[min(min(yl)) max(max(yl))])
                    end

                    print('-dtiff',['-r' num2str(imageResolution)],fullfile(saveDir, 'Text', condNames{indEffect(nInd1)}, ['All RM by ' condNames{indEffect(nInd2)}]))
                    close
                    clear ax yl dataMeans dataSD

                end
            end
        end
    end

    % 2 IND 1RM
    if numel(indEffect)>1 & numel(rmEffect)>0
        % Lines
        if statsLines
            for nRm=1:numel(rmEffect)
                for nModRm=1:numel(modalitiesRM{nRm})
                    for nInd1=1:numel(indEffect)
                        for nInd2=1:numel(indEffect)
                            if nInd1~=nInd2

                                f=figure('units','centimeters','position',[0 0 6+4*numel(modalitiesInd{nInd1}) 4+9/16*4*numel(modalitiesInd{nInd1})*numel(modalitiesInd{nInd2})],'visible','off');

                                for nMod2=1:numel(modalitiesInd{nInd2})

                                    subplot(numel(modalitiesInd{nInd2}),1,nMod2)

                                    for nMod1=1:numel(modalitiesInd{nInd1})
                                        dataMeans(nMod1)=nanmean(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nMod1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2})(:,col4means{nRm}(nModRm,:)),2));
                                        dataSD(nMod1)=nanstd(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nMod1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2})(:,col4means{nRm}(nModRm,:)),2));
                                    end


                                    for x=1:numel(modalitiesInd{nInd1})

                                        h=bar(x,dataMeans(x)); hold on
                                        h(1).FaceColor=colors{indEffect(nInd1)}(x,:);

                                        if plotSD==1
                                            if dataMeans(x)>=0
                                                SD(1)=0;
                                                SD(2)=dataSD(x);
                                            else
                                                SD(1)=dataSD(x);
                                                SD(2)=0;
                                            end
                                            errorbar(x,dataMeans(x),SD(1),SD(2),'k','LineStyle','none')
                                        end

                                        if abs(dataMeans(x))<1
                                            text(x,0.5*(dataMeans(x)),sprintf('%.3f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                        elseif abs(dataMeans(x))<10
                                            text(x,0.5*(dataMeans(x)),sprintf('%.2f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                        elseif abs(dataMeans(x))<100
                                            text(x,0.5*(dataMeans(x)),sprintf('%.1f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                        else
                                            text(x,0.5*(dataMeans(x)),sprintf('%.0f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                        end

                                    end
                                    xticks(1:numel(modalitiesInd{nInd1}))
                                    xticklabels(modalitiesInd{nInd1})
                                    if nMod2==numel(modalitiesInd{nInd2})
                                        xlabel(condNames{indEffect(nInd1)})
                                    end
                                    ylabel(units)
                                    yl(:,nMod2)=ylim;
                                    box off
                                    ax{nMod2}=gca;
                                    ax{nMod2}.XGrid='off';
                                    ax{nMod2}.YGrid='on';
                                    title(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2})

                                end

                                nAovInt3=findcolExact(aov.Effect,[condNames{indEffect(nInd1)} ':' condNames{indEffect(nInd2)} ':' condNames{rmEffect(nRm)}]);
                                if isempty(nAovInt3)
                                    nAovInt3=findComboStrings(aov.Effect,{condNames{indEffect(nInd1)}, condNames{indEffect(nInd2)}, condNames{rmEffect(nRm)}});
                                end
                                nAovInt2=findcolExact(aov.Effect,[condNames{indEffect(nInd1)} ':' condNames{indEffect(nInd2)}]);
                                if isempty(nAovInt2)
                                    nAovInt2=findcolExact(aov.Effect,[condNames{indEffect(nInd2)} ':' condNames{indEffect(nInd1)}]);
                                end
                                nAovInt=findcolExact(aov.Effect,[condNames{indEffect(nInd1)} ':' condNames{rmEffect(nRm)}]);
                                if isempty(nAovInt)
                                    nAovInt=findcolExact(aov.Effect,[condNames{indEffect(nInd1)} ':' condNames{rmEffect(nRm)}]);
                                end

                                nAov=findcolExact(aov.Effect, verifFieldName(condNames{indEffect(nInd1)}));
                                pINT3=aov{nAovInt3,6};
                                pINT2=aov{nAovInt2,6};
                                pINT=aov{nAovInt,6};
                                pMAIN=aov{nAov,6};

                                amp=[max(max(yl))-min(min(yl))];

                                for nMod2=1:numel(modalitiesInd{nInd2})

                                    for nMod1=1:numel(modalitiesInd{nInd1})
                                        dataMeans(nMod1)=nanmean(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nMod1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2})(:,col4means{nRm}(nModRm,:)),2));
                                        dataSD(nMod1)=nanstd(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nMod1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2})(:,col4means{nRm}(nModRm,:)),2));
                                    end

                                    if plotSD==1
                                        dataMeans=dataMeans+sign(dataMeans).*dataSD;
                                    end

                                    isSignificant=0;
                                    set(f,'CurrentAxes',ax{nMod2});
                                    pValues=ones(1,size(postHoc.(condNamesVerif{indEffect(nInd1)}),1));
                                    pSelected=0;
                                    if pINT3<pcritical(1) & pSelected==0
                                        isSignificant=1;
                                        phCut1=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)} 'By' condNames{rmEffect(nRm)}]){:,1}, allMod.(condNamesVerif{rmEffect(nRm)}){nModRm});
                                        phCut2=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)} 'By' condNames{rmEffect(nRm)}]){:,2}, allMod.(condNamesVerif{indEffect(nInd2)}){nMod2});
                                        phCut=intersect(phCut1, phCut2);
                                        phCut=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)} 'By' condNames{rmEffect(nRm)}])(phCut,1:end);
                                        pValues=phCut{:,7};
                                        pSelected=1;
                                    end
                                    %                                 if any([pINT<pcritical(1) pINT2<pcritical(1)]) &  pSelected==0
                                    %                                     isSignificant=1;
                                    %                                     pSelected=1;
                                    %                                     if pINT<pINT2
                                    %                                         phCut=findcolExact(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{rmEffect(nRm)}]){:,1}, allMod.(condNamesVerif{rmEffect(nRm)}){nModRm});
                                    %                                         pValues=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{rmEffect(nRm)}]){phCut,6};
                                    %                                     else
                                    %                                         phCut=findcolExact(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,1}, allMod.(condNamesVerif{indEffect(nInd2)}){nMod2});
                                    %                                         pValues=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){phCut,6};
                                    %                                     end
                                    %                                 end
                                    %                                 if pMAIN<pcritical(1) &  pSelected==0
                                    %                                     isSignificant=1;
                                    %                                     pValues=postHoc.(condNamesVerif{indEffect(nInd1)}){:,5};
                                    %                                     pSelected=1;
                                    %                                 end
                                    if pINT3<pcritical(end) & pSelected==0
                                        isSignificant=1;
                                        phCut1=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)} 'By' condNames{rmEffect(nRm)}]){:,1}, allMod.(condNamesVerif{rmEffect(nRm)}){nModRm});
                                        phCut2=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)} 'By' condNames{rmEffect(nRm)}]){:,2}, allMod.(condNamesVerif{indEffect(nInd2)}){nMod2});
                                        phCut=intersect(phCut1, phCut2);
                                        phCut=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)} 'By' condNames{rmEffect(nRm)}])(phCut,1:end);
                                        pValues=phCut{:,7};
                                        pSelected=1;
                                    end
                                    %                                 if any([pINT<pcritical(end) pINT2<pcritical(end)]) &  pSelected==0
                                    %                                     isSignificant=1;
                                    %                                     pSelected=1;
                                    %                                     if pINT<pINT2
                                    %                                         phCut=findcolExact(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{rmEffect(nRm)}]){:,1}, allMod.(condNamesVerif{rmEffect(nRm)}){nModRm});
                                    %                                         pValues=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{rmEffect(nRm)}]){phCut,6};
                                    %                                     else
                                    %                                         phCut=findcolExact(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,1}, allMod.(condNamesVerif{indEffect(nInd2)}){nMod2});
                                    %                                         pValues=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){phCut,6};
                                    %                                     end
                                    %                                 end
                                    %                                 if pMAIN<pcritical(end) &  pSelected==0
                                    %                                     isSignificant=1;
                                    %                                     pValues=postHoc.(condNamesVerif{indEffect(nInd1)}){:,5};
                                    %                                     pSelected=1;
                                    %                                 end

                                    if pSelected==1
                                        nSignificant=1;
                                        for i=1:numel(pValues)
                                            pV=pValues(i);
                                            if pV<=pcritical(end)
                                                if yl(2,nMod2)>0
                                                    if pV<pcritical(1) & isSignificant==1
                                                        hline(yl(2,nMod2)+0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNamesVerif{indEffect(nInd1)})(i,:),'lineWidth',1.5);
                                                    else
                                                        hline(yl(2,nMod2)+0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNamesVerif{indEffect(nInd1)})(i,:),'lineWidth',1.5);
                                                    end
                                                else
                                                    if pV<pcritical(1) & isSignificant==1
                                                        hline(yl(1,nMod2)-0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNamesVerif{indEffect(nInd1)})(i,:),'lineWidth',1.5);
                                                    else
                                                        hline(yl(1,nMod2)-0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNamesVerif{indEffect(nInd1)})(i,:),'lineWidth',1.5);
                                                    end
                                                end
                                                nSignificant=nSignificant+1;
                                            end
                                        end
                                    end

                                    yl(:,nMod2)=ylim;

                                end

                                for nMod2=1:numel(ax)
                                    set(ax{nMod2},'ylim',1.05*[min(min(yl)) max(max(yl))])
                                end

                                print('-dtiff',['-r' num2str(imageResolution)],fullfile(saveDir, 'Text', condNames{indEffect(nInd1)}, [condNames{rmEffect(nRm)} ' = ' allModalities{rmEffect(nRm)}{nModRm} ' By ' condNames{indEffect(nInd2)}]))
                                close
                                clear ax yl dataMeans dataSD
                            end
                        end
                    end
                end
            end
        end

        % Text
        for nRm=1:numel(rmEffect)
            for nModRm=1:numel(modalitiesRM{nRm})
                for nInd1=1:numel(indEffect)
                    for nInd2=1:numel(indEffect)
                        if nInd1~=nInd2

                            f=figure('units','centimeters','position',[0 0 6+4*numel(modalitiesInd{nInd1}) 4+9/16*4*numel(modalitiesInd{nInd1})*numel(modalitiesInd{nInd2})],'visible','off');

                            for nMod2=1:numel(modalitiesInd{nInd2})

                                subplot(numel(modalitiesInd{nInd2}),1,nMod2)

                                for nMod1=1:numel(modalitiesInd{nInd1})
                                    dataMeans(nMod1)=nanmean(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nMod1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2})(:,col4means{nRm}(nModRm,:)),2));
                                    dataSD(nMod1)=nanstd(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nMod1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2})(:,col4means{nRm}(nModRm,:)),2));
                                end

                                for x=1:numel(modalitiesInd{nInd1})

                                    h=bar(x,dataMeans(x)); hold on
                                    h(1).FaceColor=colors{indEffect(nInd1)}(x,:);

                                    if plotSD==1
                                        if dataMeans(x)>=0
                                            SD(1)=0;
                                            SD(2)=dataSD(x);
                                        else
                                            SD(1)=dataSD(x);
                                            SD(2)=0;
                                        end
                                        errorbar(x,dataMeans(x),SD(1),SD(2),'k','LineStyle','none')
                                    end

                                    if abs(dataMeans(x))<1
                                        text(x,0.5*(dataMeans(x)),sprintf('%.3f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                    elseif abs(dataMeans(x))<10
                                        text(x,0.5*(dataMeans(x)),sprintf('%.2f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                    elseif abs(dataMeans(x))<100
                                        text(x,0.5*(dataMeans(x)),sprintf('%.1f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                    else
                                        text(x,0.5*(dataMeans(x)),sprintf('%.0f',dataMeans(x)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color',rgb('white'),'FontWeight','bold')
                                    end

                                end
                                xticks(1:numel(modalitiesInd{nInd1}))
                                xticklabels(modalitiesInd{nInd1})
                                if nMod2==numel(modalitiesInd{nInd2})
                                    xlabel(condNames{indEffect(nInd1)})
                                end
                                ylabel(units)
                                yl(:,nMod2)=ylim;
                                box off
                                ax{nMod2}=gca;
                                ax{nMod2}.XGrid='off';
                                ax{nMod2}.YGrid='on';
                                title(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2})

                            end

                            nAovInt3=findcolExact(aov.Effect,[condNames{indEffect(nInd1)} ':' condNames{indEffect(nInd2)} ':' condNames{rmEffect(nRm)}]);
                            if isempty(nAovInt3)
                                nAovInt3=findComboStrings(aov.Effect,{condNames{indEffect(nInd1)}, condNames{indEffect(nInd2)}, condNames{rmEffect(nRm)}});
                            end
                            nAovInt2=findcolExact(aov.Effect,[condNames{indEffect(nInd1)} ':' condNames{indEffect(nInd2)}]);
                            if isempty(nAovInt2)
                                nAovInt2=findcolExact(aov.Effect,[condNames{indEffect(nInd2)} ':' condNames{indEffect(nInd1)}]);
                            end
                            nAovInt=findcolExact(aov.Effect,[condNames{indEffect(nInd1)} ':' condNames{rmEffect(nRm)}]);
                            if isempty(nAovInt)
                                nAovInt=findcolExact(aov.Effect,[condNames{indEffect(nInd1)} ':' condNames{rmEffect(nRm)}]);
                            end

                            nAov=findcolExact(aov.Effect, verifFieldName(condNames{indEffect(nInd1)}));
                            pINT3=aov{nAovInt3,6};
                            pINT2=aov{nAovInt2,6};
                            pINT=aov{nAovInt,6};
                            pMAIN=aov{nAov,6};

                            amp=[max(max(yl))-min(min(yl))];

                            for nMod2=1:numel(modalitiesInd{nInd2})

                                for nMod1=1:numel(modalitiesInd{nInd1})
                                    dataMeans(nMod1)=nanmean(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nMod1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2})(:,col4means{nRm}(nModRm,:)),2));
                                    dataSD(nMod1)=nanstd(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNamesVerif{indEffect(nInd1)}){nMod1}).(allMod.(condNamesVerif{indEffect(nInd2)}){nMod2})(:,col4means{nRm}(nModRm,:)),2));
                                end

                                dataMeans4pv=dataMeans;
                                if plotSD==1
                                    dataMeans4pv=dataMeans+sign(dataMeans).*dataSD;
                                end

                                isSignificant=0;
                                set(f,'CurrentAxes',ax{nMod2});
                                pValues=ones(1,size(postHoc.(condNamesVerif{indEffect(nInd1)}),1));
                                pSelected=0;
                                if pINT3<pcritical(1) & pSelected==0
                                    isSignificant=1;
                                    phCut1=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)} 'By' condNames{rmEffect(nRm)}]){:,1}, allMod.(condNamesVerif{rmEffect(nRm)}){nModRm});
                                    phCut2=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)} 'By' condNames{rmEffect(nRm)}]){:,2}, allMod.(condNamesVerif{indEffect(nInd2)}){nMod2});
                                    phCut=intersect(phCut1, phCut2);
                                    phCut=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)} 'By' condNames{rmEffect(nRm)}])(phCut,1:end);
                                    pValues=phCut{:,7};
                                    pSelected=1;
                                end
                                %                             if any([pINT<pcritical(1) pINT2<pcritical(1)]) &  pSelected==0
                                %                                 isSignificant=1;
                                %                                 pSelected=1;
                                %                                 if pINT<pINT2
                                %                                     phCut=findcolExact(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{rmEffect(nRm)}]){:,1}, allMod.(condNamesVerif{rmEffect(nRm)}){nModRm});
                                %                                     pValues=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{rmEffect(nRm)}]){phCut,6};
                                %                                 else
                                %                                     phCut=findcolExact(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,1}, allMod.(condNamesVerif{indEffect(nInd2)}){nMod2});
                                %                                     pValues=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){phCut,6};
                                %                                 end
                                %                             end
                                %                             if pMAIN<pcritical(1) &  pSelected==0
                                %                                 isSignificant=1;
                                %                                 pValues=postHoc.(condNamesVerif{indEffect(nInd1)}){:,5};
                                %                                 pSelected=1;
                                %                             end
                                if pINT3<pcritical(end) & pSelected==0
                                    isSignificant=1;
                                    phCut1=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)} 'By' condNames{rmEffect(nRm)}]){:,1}, allMod.(condNamesVerif{rmEffect(nRm)}){nModRm});
                                    phCut2=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)} 'By' condNames{rmEffect(nRm)}]){:,2}, allMod.(condNamesVerif{indEffect(nInd2)}){nMod2});
                                    phCut=intersect(phCut1, phCut2);
                                    phCut=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)} 'By' condNames{rmEffect(nRm)}])(phCut,1:end);
                                    pValues=phCut{:,7};
                                    pSelected=1;
                                end
                                %                             if any([pINT<pcritical(end) pINT2<pcritical(end)]) &  pSelected==0
                                %                                 isSignificant=1;
                                %                                 pSelected=1;
                                %                                 if pINT<pINT2
                                %                                     phCut=findcolExact(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{rmEffect(nRm)}]){:,1}, allMod.(condNamesVerif{rmEffect(nRm)}){nModRm});
                                %                                     pValues=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{rmEffect(nRm)}]){phCut,6};
                                %                                 else
                                %                                     phCut=findcolExact(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,1}, allMod.(condNamesVerif{indEffect(nInd2)}){nMod2});
                                %                                     pValues=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){phCut,6};
                                %                                 end
                                %                             end
                                %                             if pMAIN<pcritical(end) &  pSelected==0
                                %                                 isSignificant=1;
                                %                                 pValues=postHoc.(condNamesVerif{indEffect(nInd1)}){:,5};
                                %                                 pSelected=1;
                                %                             end

                                nColSignificant=ones(1,numel(dataMeans));
                                nSignificant=1;
                                if pSelected==1
                                    for i=1:numel(pValues)
                                        pV=pValues(i);
                                        if pV<=pcritical(end)
                                            if abs(yl(2))>abs(yl(1))
                                                if abs(dataMeans(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1)))>abs(dataMeans(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2)))
                                                    addPvalue(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1), dataMeans4pv(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1))+0.0225*numel(modalitiesInd{nInd2})*nColSignificant(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1))*amp, pV, pcritical, colors{indEffect(nInd1)}(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2),:))
                                                    nColSignificant(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1))=nColSignificant(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1))+1;
                                                else
                                                    addPvalue(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2), dataMeans4pv(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2))+0.0225*numel(modalitiesInd{nInd2})*nColSignificant(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2))*amp, pV, pcritical, colors{indEffect(nInd1)}(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1),:))
                                                    nColSignificant(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2))=nColSignificant(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2))+1;

                                                end
                                            else
                                                if abs(dataMeans(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1)))>abs(dataMeans(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2)))
                                                    addPvalue(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1), dataMeans4pv(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1))-0.0225*numel(modalitiesInd{nInd2})*nColSignificant(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1))*amp, pV, pcritical, colors{indEffect(nInd1)}(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2),:))
                                                    nColSignificant(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1))=nColSignificant(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1))+1;
                                                else
                                                    addPvalue(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2), dataMeans4pv(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2))-0.0225*numel(modalitiesInd{nInd2})*nColSignificant(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2))*amp, pV, pcritical, colors{indEffect(nInd1)}(order4ES.(condNamesVerif{indEffect(nInd1)})(i,1),:))
                                                    nColSignificant(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2))=nColSignificant(order4ES.(condNamesVerif{indEffect(nInd1)})(i,2))+1;
                                                end
                                            end
                                            nSignificant=nSignificant+1;
                                        end
                                    end
                                end
                            end

                            yl=yl*1.05.^max(nColSignificant);
                            for nMod2=1:numel(ax)
                                set(ax{nMod2},'ylim',1.05*[min(min(yl)) max(max(yl))])
                            end

                            print('-dtiff',['-r' num2str(imageResolution)],fullfile(saveDir, 'Text', condNames{indEffect(nInd1)}, [condNames{rmEffect(nRm)} ' = ' allModalities{rmEffect(nRm)}{nModRm} ' By ' condNames{indEffect(nInd2)}]))
                            close
                            clear ax yl dataMeans dataSD

                        end
                    end
                end
            end
        end
    end
end

end