% create graph and tables for a 4way anova with 2independant effect and 2rm depending on the
% project parameters values

function [Tables]=fctAOV(data,stats,display,units,saveDir)

%% Unpack structures and create save directory
unitsTmp=units;
warning('off');
mkdir(saveDir)
unpackStruct(stats)
unpackStruct(display)
units=unitsTmp;

%% Statistical names and order based on the names
% modalities
loopRm=1;
for nCond=1:numel(cond4effect)
    allModalities{nCond}=unique(cond4effect{nCond},'stable');
    allMod.(condNames{nCond})=unique(cond4effect{nCond},'stable');
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

if sum(isRM)==1
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
                order4ES.([condNames{rmEffect(nCond1)} 'By' condNames{rmEffect(nCond2)}])=order4EStmp;
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
                varData(s,j)=nanmean(data{s,i,j});
            end
        end
        t{i}=varData;
        clear varData
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create tables for means
if numel(indEffect)==1
    cNames=['ID' condNames{indEffect(1)} condNames4Table];
elseif numel(indEffect)==2
    cNames=['ID' condNames{indEffect(1)} condNames{indEffect(2)} condNames4Table];
end

tXl=[];
for i=1:numel(cond4effect{rmEffect(1)})
    tXl=[tXl t{i}];
end
% nan the subjects with issues
tXl(isnan(mean(tXl,2)),:)=nan;
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
    tXLmeans(3,numel(indEffect)+1+i)=table(nan);
end
data4plot.allData=tXl;
loop=4;
if numel(indEffect)>0
    for nInd1=1:numel(indEffect)
        for nMod1=1:numel(modalitiesInd{nInd1})
            for i=1:numel(condNames4Table)
                tXLmeans(loop,numel(indEffect)+1+i)=table(nanmean(tXl(idxIndependantEffect(:,nInd1)==nMod1,i)));
                tXLmeans(loop+1,numel(indEffect)+1+i)=table(nanstd(tXl(idxIndependantEffect(:,nInd1)==nMod1,i)));
                tXLmeans(loop+2,numel(indEffect)+1+i)=table(nan);
            end
            data4plot.(condNames{[indEffect(nInd1)]}).(allModalities{[indEffect(nInd1)]}{nMod1})=tXl(idxIndependantEffect(:,nInd1)==nMod1,:);
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
                        data4plot.([condNames{[indEffect(nInd1)]} 'By' condNames{[indEffect(nInd2)]}]).([allModalities{[indEffect(nInd1)]}{nMod1} allModalities{[indEffect(nInd2)]}{nMod2}])=tXl(idxIndependantEffect(:,nInd1)==nMod1 & idxIndependantEffect(:,nInd2)==nMod2,:);
                        data4plot.([condNames{[indEffect(nInd1)]} 'By' condNames{[indEffect(nInd2)]}]).(allModalities{[indEffect(nInd1)]}{nMod1}).(allModalities{indEffect(nInd2)}{nMod2})=tXl(idxIndependantEffect(:,nInd1)==nMod1 & idxIndependantEffect(:,nInd2)==nMod2,:);
                        data4plot.([condNames{[indEffect(nInd2)]} 'By' condNames{[indEffect(nInd1)]}]).(allModalities{[indEffect(nInd2)]}{nMod2}).(allModalities{indEffect(nInd1)}{nMod1})=tXl(idxIndependantEffect(:,nInd1)==nMod1 & idxIndependantEffect(:,nInd2)==nMod2,:);
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
                                    data4plot.([condNames{[indEffect(nInd1)]} 'By' condNames{[indEffect(nInd2)]} 'By' condNames{[indEffect(nInd3)]}]).([allModalities{[indEffect(nInd1)]}{nMod1} allModalities{[indEffect(nInd2)]}{nMod2} allModalities{[indEffect(nInd3)]}{nMod3}])=tXl(idxIndependantEffect(:,nInd1)==nMod1 & idxIndependantEffect(:,nInd2)==nMod2 & idxIndependantEffect(:,nInd3)==nMod3,:);
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
for k=1:numel(cond4effect{rmEffect(1)})
    data4stats(:,k,:)=t{k};
end

if numel(rmEffect)==2 & numel(indEffect)==2

    [tbl,rm]=mixed_anova(data4stats,idxIndependantEffect, {condNames{rmEffect(1)}, condNames{rmEffect(2)}}, {condNames{indEffect(1)}, condNames{indEffect(2)}});
    rm.WithinDesign.(condNames{rmEffect(1)})=effectRMaov{1};
    rm.WithinDesign.(condNames{rmEffect(2)})=effectRMaov{2};

elseif numel(rmEffect)==1 & numel(indEffect)==1

    [tbl,rm]=mixed_anova(data4stats,idxIndependantEffect, {condNames{rmEffect(1)}}, {condNames{indEffect(1)}});
    rm.WithinDesign.(condNames{rmEffect(1)})=effectRMaov{1};

elseif numel(rmEffect)==2 & numel(indEffect)==1

    [tbl,rm]=mixed_anova(data4stats,idxIndependantEffect, {condNames{rmEffect(1)}, condNames{rmEffect(2)}}, {condNames{indEffect(1)}});
    rm.WithinDesign.(condNames{rmEffect(1)})=effectRMaov{1};
    rm.WithinDesign.(condNames{rmEffect(2)})=effectRMaov{2};

end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1 effect
for nCond=1:numel(condNames)
    ph=multcompare(rm,condNames{nCond},'ComparisonType',postHocType);
    [whichPH,order4ES.(condNames{nCond})]=findPH(allModalities{nCond});
    postHoc.(condNames{nCond})=ph(whichPH,:);
    if any(nCond==indEffect)
        postHoc.(condNames{nCond}).([condNames{nCond} '_1'])=modalitiesInd{nCond}(postHoc.(condNames{nCond}).([condNames{nCond} '_1']));
        postHoc.(condNames{nCond}).([condNames{nCond} '_2'])=modalitiesInd{nCond}(postHoc.(condNames{nCond}).([condNames{nCond} '_2']));
    end
end

%% 2 effects
if numel(condNames)>1
    for nCond1=1:numel(condNames)
        for nCond2=1:numel(condNames)
            if nCond2~=nCond1
                ph=multcompare(rm,condNames{nCond1},'By',condNames{nCond2},'ComparisonType',postHocType);
                [whichPH]=findPHint(allModalities{nCond1},allModalities{nCond2});
                postHoc.([condNames{nCond1} 'By'  condNames{nCond2}])=ph(whichPH,:);
                if any(nCond1==indEffect)
                    postHoc.([condNames{nCond1} 'By'  condNames{nCond2}]).([condNames{nCond1} '_1'])=modalitiesInd{nCond1}(postHoc.([condNames{nCond1} 'By'  condNames{nCond2}]).([condNames{nCond1} '_1']));
                    postHoc.([condNames{nCond1} 'By'  condNames{nCond2}]).([condNames{nCond1} '_2'])=modalitiesInd{nCond1}(postHoc.([condNames{nCond1} 'By'  condNames{nCond2}]).([condNames{nCond1} '_2']));
                end
                if any(nCond2==indEffect)
                    postHoc.([condNames{nCond1} 'By'  condNames{nCond2}]).([condNames{nCond2}])=modalitiesInd{nCond2}(postHoc.([condNames{nCond1} 'By'  condNames{nCond2}]).([condNames{nCond2}]));
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
                        [tbl,rm]=mixed_anova(tXl,idxIndependantEffect, {[condNames{nCond2} 'By' condNames{nCond3}]}, {condNames{indEffect(1)}, condNames{indEffect(2)}});
                    elseif numel(condNames)==3 & any(nCond1==indEffect)
                        [tbl,rm]=mixed_anova(tXl,idxIndependantEffect, {[condNames{nCond2} 'By' condNames{nCond3}]}, {condNames{indEffect(nCond1)}});
                    end

                    rm.WithinDesign.([condNames{nCond2} 'By' condNames{nCond3}])=transpose(string(condNames4Table));
                    ph=multcompare(rm,condNames{nCond1},'By', [condNames{nCond2} 'By' condNames{nCond3}], 'ComparisonType',postHocType);
                    [whichPH]=findPHint(allModalities{nCond1},condNames4Table);
                    if any(nCond1==indEffect)
                        ph.([condNames{nCond1} '_1'])=modalitiesInd{nCond1}(ph.([condNames{nCond1} '_1']));
                        ph.([condNames{nCond1} '_2'])=modalitiesInd{nCond1}(ph.([condNames{nCond1} '_2']));
                    end
                    ph=ph(whichPH,:);
                    phNew=formatPH3_firstCol(ph,allModalities([nCond2 nCond3]), condNames([nCond2 nCond3]));
                    postHoc.([condNames{nCond1} 'By'  condNames{nCond2} 'By'  condNames{nCond3}])=phNew;

                    % ind*ind*ind/rm
                elseif all([nCond3~=nCond2 & nCond3~=nCond1 & nCond2~=nCond1 any(nCond1==indEffect) any(nCond2==indEffect)])

                    for s=1:size(idxIndependantEffect,1)
                        int2.([condNames{nCond1} 'By' condNames{nCond2}])(s,1)=string([cond4effect{nCond1}{s} '' cond4effect{nCond2}{s}]);
                    end
                    modInt2.([condNames{nCond1} 'By' condNames{nCond2}])=unique(int2.([condNames{nCond1} 'By' condNames{nCond2}]));
                    allMod.([condNames{nCond1} 'By' condNames{nCond2}])=unique(cellstr(int2.([condNames{nCond1} 'By' condNames{nCond2}])));
                    for s=1:numel(int2.([condNames{nCond1} 'By' condNames{nCond2}]))
                        indInt2.([condNames{nCond1} 'By' condNames{nCond2}])(s,1)=findcol(modInt2.([condNames{nCond1} 'By' condNames{nCond2}]), int2.([condNames{nCond1} 'By' condNames{nCond2}])(s));
                    end
                    [tbl,rm]=mixed_anova(data4stats,indInt2.([condNames{nCond1} 'By' condNames{nCond2}]), {condNames{rmEffect(1)}, condNames{rmEffect(2)}}, {[condNames{nCond1} 'By' condNames{nCond2}]});
                    rm.WithinDesign.(condNames{rmEffect(1)})=effectRMaov{1};
                    rm.WithinDesign.(condNames{rmEffect(2)})=effectRMaov{2};
                    ph=multcompare(rm,[condNames{nCond1} 'By' condNames{nCond2}], 'By', condNames{nCond3}, 'ComparisonType',postHocType);
                    [whichPH]=findPHint(modInt2.([condNames{nCond1} 'By' condNames{nCond2}]),allModalities{nCond3});
                    [~, order4ES.([condNames{nCond1} 'By' condNames{nCond2}])]=findPH(modInt2.([condNames{nCond1} 'By' condNames{nCond2}]));
                    ph=ph(whichPH,:);
                    ph.([condNames{nCond1} 'By' condNames{nCond2} '_1'])=modInt2.([condNames{nCond1} 'By' condNames{nCond2}])(ph.([condNames{nCond1} 'By' condNames{nCond2} '_1']));
                    ph.([condNames{nCond1} 'By' condNames{nCond2} '_2'])=modInt2.([condNames{nCond1} 'By' condNames{nCond2}])(ph.([condNames{nCond1} 'By' condNames{nCond2} '_2']));
                    phNew=formatPH3v2(ph,allModalities([nCond1 nCond2 nCond3]), condNames([nCond1 nCond2 nCond3]));
                    postHoc.([condNames{nCond1} 'By'  condNames{nCond2} 'By'  condNames{nCond3}])=phNew;

                    %ind/rm*ind*ind
                elseif all([nCond3>nCond2 & nCond3~=nCond1 & nCond2~=nCond1 any(nCond2==indEffect) any(nCond3==indEffect)])

                    for s=1:size(idxIndependantEffect,1)
                        int2.([condNames{nCond2} 'By' condNames{nCond3}])(s,1)=string([cond4effect{nCond2}{s} '' cond4effect{nCond3}{s}]);
                    end
                    modInt2.([condNames{nCond2} 'By' condNames{nCond3}])=unique(int2.([condNames{nCond2} 'By' condNames{nCond3}]));
                    for s=1:numel( int2.([condNames{nCond2} 'By' condNames{nCond3}]))
                        indInt2.([condNames{nCond2} 'By' condNames{nCond3}])(s,1)=findcol(modInt2.([condNames{nCond2} 'By' condNames{nCond3}]), int2.([condNames{nCond2} 'By' condNames{nCond3}])(s));
                    end
                    [tbl,rm]=mixed_anova(data4stats,indInt2.([condNames{nCond2} 'By' condNames{nCond3}]), {condNames{rmEffect(1)}, condNames{rmEffect(2)}}, {[condNames{nCond2} 'By' condNames{nCond3}]});
                    rm.WithinDesign.(condNames{rmEffect(1)})=effectRMaov{1};
                    rm.WithinDesign.(condNames{rmEffect(2)})=effectRMaov{2};
                    ph=multcompare(rm,condNames{nCond1},'By', [condNames{nCond2} 'By' condNames{nCond3}], 'ComparisonType',postHocType);
                    [whichPH]=findPHint(allModalities{nCond1},modInt2.([condNames{nCond2} 'By' condNames{nCond3}]));
                    ph=ph(whichPH,:);
                    ph.([condNames{nCond2} 'By' condNames{nCond3}])=modInt2.([condNames{nCond2} 'By' condNames{nCond3}])(ph.([condNames{nCond2} 'By' condNames{nCond3}]));
                    phNew=formatPH3_firstCol(ph,allModalities([nCond2 nCond3]), condNames([nCond2 nCond3]));
                    postHoc.([condNames{nCond1} 'By'  condNames{nCond2} 'By'  condNames{nCond3}])=phNew;

                    % rm*rm*ind
                elseif all([nCond3~=nCond2 & nCond3~=nCond1 & nCond2~=nCond1 any(nCond1==rmEffect) any(nCond2==rmEffect)]) & any(nCond3==indEffect)

                    [tbl,rm]=mixed_anova(tXl,idxIndependantEffect(:,indEffect(nCond3)), {[condNames{nCond1} 'By' condNames{nCond2}]}, {condNames{indEffect(nCond3)}});
                    rm.WithinDesign.([condNames{nCond1} 'By' condNames{nCond2}])=condNames4Table';
                    ph=multcompare(rm, [condNames{nCond1} 'By' condNames{nCond2}], 'By', condNames{nCond3}, 'ComparisonType',postHocType);
                    phNew=formatPH3(ph,allModalities([nCond1 nCond2 nCond3]), condNames([nCond1 nCond2 nCond3]));
                    postHoc.([condNames{nCond1} 'By'  condNames{nCond2} 'By'  condNames{nCond3}])=phNew;

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
                            int2.([condNames{nCond1} 'By' condNames{nCond2}])(s,1)=string([cond4effect{nCond1}{s} ' & ' cond4effect{nCond2}{s}]);
                        end
                        modInt2.([condNames{nCond1} 'By' condNames{nCond2}])=unique(int2.([condNames{nCond1} 'By' condNames{nCond2}]));
                        for s=1:numel(int2.([condNames{nCond1} 'By' condNames{nCond2}]))
                            indInt2.([condNames{nCond1} 'By' condNames{nCond2}])(s,1)=findcol(modInt2.([condNames{nCond1} 'By' condNames{nCond2}]),int2.([condNames{nCond1} 'By' condNames{nCond2}])(s));
                        end
                        [tbl,rm]=mixed_anova(tXl,indInt2.([condNames{nCond1} 'By' condNames{nCond2}]), {[condNames{nCond3} 'By' condNames{nCond4}]}, {[condNames{nCond1} 'By' condNames{nCond2}]});
                        rm.WithinDesign.([condNames{nCond3} 'By' condNames{nCond4}])=transpose(string(condNames4Table));
                        ph=multcompare(rm,[condNames{nCond1} 'By' condNames{nCond2}],'By', [condNames{nCond3} 'By' condNames{nCond4}], 'ComparisonType',postHocType);
                        [whichPH]=findPHint(modInt2.([condNames{nCond1} 'By' condNames{nCond2}]),condNames4Table);
                        phName=[condNames{nCond1}  'By'  condNames{nCond2} 'By' condNames{nCond3} 'By'  condNames{nCond4}];
                        ind2replaceName=[condNames{nCond1} 'By' condNames{nCond2}];
                        ph=ph(whichPH,:);
                        ph.([ind2replaceName '_1'])=modInt2.([condNames{nCond1} 'By' condNames{nCond2}])(ph.([ind2replaceName '_1']));
                        ph.([ind2replaceName '_2'])=modInt2.([condNames{nCond1} 'By' condNames{nCond2}])(ph.([ind2replaceName '_2']));
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
                            int2.([condNames{nCond1} 'By' condNames{nCond2}])(s,1)=string([cond4effect{nCond1}{s} ' & ' cond4effect{nCond2}{s}]);
                        end
                        modInt2.([condNames{nCond1} 'By' condNames{nCond2}])=unique(int2.([condNames{nCond1} 'By' condNames{nCond2}]));
                        for s=1:numel(int2.([condNames{nCond1} 'By' condNames{nCond2}]))
                            indInt2.([condNames{nCond1} 'By' condNames{nCond2}])(s,1)=findcol(modInt2.([condNames{nCond1} 'By' condNames{nCond2}]),int2.([condNames{nCond1} 'By' condNames{nCond2}])(s));
                        end

                        [tbl,rm]=mixed_anova(tXl,indInt2.([condNames{nCond1} 'By' condNames{nCond2}]), {[condNames{nCond3} 'By' condNames{nCond4}]}, {[condNames{nCond1} 'By' condNames{nCond2}]});
                        rm.WithinDesign.([condNames{nCond3} 'By' condNames{nCond4}])=transpose(string(condNames4Table));

                        ph=multcompare(rm,[condNames{nCond3} 'By' condNames{nCond4}],'By', [condNames{nCond1} 'By' condNames{nCond2}], 'ComparisonType',postHocType);
                        [whichPH]=findPHint(condNames4Table, modInt2.([condNames{nCond1} 'By' condNames{nCond2}]));
                        phName=[condNames{nCond3}  'By'  condNames{nCond4} 'By' condNames{nCond1} 'By'  condNames{nCond2}];
                        ind2replaceName=[condNames{nCond1} 'By' condNames{nCond2}];
                        ph=ph(whichPH,:);
                        ph.(ind2replaceName)=modInt2.([condNames{nCond1} 'By' condNames{nCond2}])(ph.(ind2replaceName));

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

    tablemeans.(condNames{rmEffect(nCond)})=[tempTable; tempTableMeans];
    if  numel(indEffect)==1
        tablemeans.(condNames{rmEffect(nCond)}).Properties.VariableNames=['var' condNames{indEffect(1)} cond4effect{rmEffect(nCond)}];
    elseif numel(indEffect)==2
        tablemeans.(condNames{rmEffect(nCond)}).Properties.VariableNames=['var' condNames{indEffect(1)} condNames{indEffect(2)} cond4effect{rmEffect(nCond)}];
    elseif numel(indEffect)==3
        tablemeans.(condNames{rmEffect(nCond)}).Properties.VariableNames=['var' condNames{indEffect(1)} condNames{indEffect(2)} condNames{indEffect(3)} cond4effect{rmEffect(nCond)}];
    end
    clear tempTableMeans tempTable


    % Effect sizes
    phFieldnames=fieldnames(postHoc);
    % simple rm effects
    for o=1:size(order4ES.(condNames{rmEffect(nCond)}),1)
        dES{1}=means(:,order4ES.(condNames{rmEffect(nCond)})(o,1));
        dES{2}=means(:,order4ES.(condNames{rmEffect(nCond)})(o,2));
        ES(o,1)=esCalculation0D(dES);
    end
    postHoc.(condNames{rmEffect(nCond)})=[postHoc.(condNames{rmEffect(nCond)}) table(ES)];
    clear ES

    % interaction of rm with independant effect
    if numel(indEffect)>0
        for nInd=1:numel(indEffect)
            [~,b]=sort(allModalities{indEffect(nInd)});
            for nMod=1:numel(modalitiesInd{nInd})
                for o=1:size(order4ES.(condNames{rmEffect(nCond)}),1)
                    dES{1}=means(idxIndependantEffect(:,nInd)==nMod,order4ES.(condNames{rmEffect(nCond)})(o,1));
                    dES{2}=means(idxIndependantEffect(:,nInd)==nMod,order4ES.(condNames{rmEffect(nCond)})(o,2));
                    ES(o,b(nMod))=esCalculation0D(dES);
                end
            end
            ES=ES(:);
            postHoc.([condNames{rmEffect(nCond)} 'By' condNames{indEffect(nInd)}])=[postHoc.([condNames{rmEffect(nCond)} 'By' condNames{indEffect(nInd)}]) table(ES)];
            clear ES
        end
    end

    % interaction of rm with douple independant
    if numel(indEffect)>1
        fieldNames=fieldnames(modInt2);
        for nInd=1:numel(fieldNames)
            [~,b]=sort(modInt2.(fieldNames{nInd}));
            for nMod=1:max(b)
                for o=1:size(order4ES.(condNames{rmEffect(nCond)}),1)
                    dES{1}=means(indInt2.(fieldNames{nInd})==nMod,order4ES.(condNames{rmEffect(nCond)})(o,1));
                    dES{2}=means(indInt2.(fieldNames{nInd})==nMod,order4ES.(condNames{rmEffect(nCond)})(o,2));
                    ES(o,b(nMod))=esCalculation0D(dES);
                end
            end
            ES=ES(:);
            if ~isempty(findcolExact(phFieldnames,[condNames{rmEffect(nCond)} 'By' fieldNames{nInd}]))
                postHoc.([condNames{rmEffect(nCond)} 'By' fieldNames{nInd}])=[postHoc.([condNames{rmEffect(nCond)} 'By' fieldNames{nInd}]) table(ES)];
            end
            clear ES
        end
    end
end

% interaction betwween rm effects
if numel(rmEffect)>1
    for nCond1=1:sum(isRM)
        for nCond2=1:sum(isRM)
            if nCond1~=nCond2
                for o=1:size(order4ES.([condNames{rmEffect(nCond1)} 'By'  condNames{rmEffect(nCond2)}]),1)
                    dES{1}=tXl(:,order4ES.([condNames{rmEffect(nCond1)} 'By'  condNames{rmEffect(nCond2)}])(o,1));
                    dES{2}=tXl(:,order4ES.([condNames{rmEffect(nCond1)} 'By'  condNames{rmEffect(nCond2)}])(o,2));
                    ES(o,1)=esCalculation0D(dES);
                end
                postHoc.([condNames{rmEffect(nCond1)} 'By'  condNames{rmEffect(nCond2)}])=[postHoc.([condNames{rmEffect(nCond1)} 'By'  condNames{rmEffect(nCond2)}]) table(ES)];
                clear ES
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
                        for o=1:size(order4ES.([condNames{rmEffect(nCond1)} 'By'  condNames{rmEffect(nCond2)}]),1)
                            dES{1}=tXl(idxIndependantEffect(:,nInd)==nModInd,order4ES.([condNames{rmEffect(nCond1)} 'By'  condNames{rmEffect(nCond2)}])(o,1));
                            dES{2}=tXl(idxIndependantEffect(:,nInd)==nModInd,order4ES.([condNames{rmEffect(nCond1)} 'By'  condNames{rmEffect(nCond2)}])(o,2));
                            ES(o,nModInd)=esCalculation0D(dES);
                        end
                    end
                    ES=ES(:);
                    postHoc.([condNames{rmEffect(nCond1)} 'By'  condNames{rmEffect(nCond2)} 'By' condNames{indEffect(nInd)}])=[postHoc.([condNames{rmEffect(nCond1)} 'By'  condNames{rmEffect(nCond2)} 'By' condNames{indEffect(nInd)}]) table(ES)];
                    clear ES
                end
            end
        end
    end
end


% simple independant effect
dataAllRm=nanmean(tXl,2);
if numel(indEffect)>0
    for nInd=1:numel(indEffect)
        for o=1:size(order4ES.(condNames{indEffect(nInd)}),1)
            dES{1}=dataAllRm(order4ES.(condNames{indEffect(nInd)})(o,1)==idxIndependantEffect(:,nInd));
            dES{2}=dataAllRm(order4ES.(condNames{indEffect(nInd)})(o,2)==idxIndependantEffect(:,nInd));
            ES(o,1)=esCalculation0D(dES);
        end
        postHoc.(condNames{indEffect(nInd)})=[postHoc.(condNames{indEffect(nInd)}) table(ES)];
        clear ES
    end
end

% interaction of indpendant with rm effect
if numel(indEffect)>0 & numel(rmEffect)>0
    for nInd=1:numel(indEffect)
        for nRm=1:numel(rmEffect)
            for o=1:size(order4ES.(condNames{indEffect(nInd)}),1)
                for oRm=1:size(col4means{nRm},1)
                    dES{1}=nanmean(tXl(order4ES.(condNames{indEffect(nInd)})(o,1)==idxIndependantEffect(:,nInd),col4means{nRm}(oRm,:)),2);
                    dES{2}=nanmean(tXl(order4ES.(condNames{indEffect(nInd)})(o,2)==idxIndependantEffect(:,nInd),col4means{nRm}(oRm,:)),2);
                    ES(o,oRm)=esCalculation0D(dES);
                end
            end
            ES=ES(:);
            postHoc.([condNames{indEffect(nInd)} 'By' condNames{rmEffect(nRm)}])=[postHoc.([condNames{indEffect(nInd)} 'By' condNames{rmEffect(nRm)}]) table(ES)];
            clear ES
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
                    for o=1:size(order4ES.(condNames{indEffect(nInd)}),1)
                        for oRm=1:size(tXl,2)
                            dES{1}=(tXl(order4ES.(condNames{indEffect(nInd)})(o,1)==idxIndependantEffect(:,nInd),oRm));
                            dES{2}=(tXl(order4ES.(condNames{indEffect(nInd)})(o,2)==idxIndependantEffect(:,nInd),oRm));
                            ES(o,oRm)=esCalculation0D(dES);
                        end
                    end
                    ES=ES(:);
                    postHoc.([condNames{indEffect(nInd)} 'By' condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}])=[postHoc.([condNames{indEffect(nInd)} 'By' condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}]) table(ES)];
                    clear ES
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
                    for o=1:size(order4ES.(condNames{indEffect(nInd1)}),1)
                        dES{1}=data4ES(order4ES.(condNames{indEffect(nInd1)})(o,1)==idxIndependantEffect(:,nInd1) & idxIndependantEffect(:,nInd2)==nMod2);
                        dES{2}=data4ES(order4ES.(condNames{indEffect(nInd1)})(o,2)==idxIndependantEffect(:,nInd1) & idxIndependantEffect(:,nInd2)==nMod2);
                        ES(o,b(nMod2))=esCalculation0D(dES);
                    end
                end
                ES=ES(:);
                postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}])=[postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]) table(ES)];
                clear ES
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
                            for o=1:size(order4ES.(condNames{nInd1}),1)
                                dES{1}=data4ES(order4ES.(condNames{indEffect(nInd1)})(o,1)==idxIndependantEffect(:,nInd1) & idxIndependantEffect(:,nInd2)==nMod2);
                                dES{2}=data4ES(order4ES.(condNames{indEffect(nInd1)})(o,2)==idxIndependantEffect(:,nInd1) & idxIndependantEffect(:,nInd2)==nMod2);
                                ES(o,nMod2,bRm(nModRm))=esCalculation0D(dES);
                            end
                        end
                    end
                    ES=ES(:);
                    if ~isempty(findcolExact(phFieldnames,[condNames{nInd1} 'By' condNames{nInd2} 'By' condNames{rmEffect(nRm)}]))
                        postHoc.([condNames{nInd1} 'By' condNames{nInd2} 'By' condNames{rmEffect(nRm)}])=[postHoc.([condNames{nInd1} 'By' condNames{nInd2} 'By' condNames{rmEffect(nRm)}]) table(ES)];
                    end
                    clear ES
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
writetable(aov,fullfile(saveDir, 'Tables.xlsx'),"Sheet",'AOV')
for i=1:numel(fieldNames4SaveMeans)
    writetable(tablemeans.(fieldNames4SaveMeans{i}),fullfile(saveDir, 'Tables.xlsx'),"Sheet",['MEANS_' fieldNames4SaveMeans{i}])
end
for i=1:numel(fieldNames4SavePH)
    writetable(postHoc.(fieldNames4SavePH{i}),fullfile(saveDir, 'Tables.xlsx'),"Sheet",['PH_' fieldNames4SavePH{i}])
end
Tables.means=tablemeans;
Tables.aov=aov;
Tables.posthoc=postHoc;
save(fullfile(saveDir, 'Tables'),'Tables')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PLOT
for nCond=1:numel(condNames)
    if plotLines
        mkdir(fullfile(saveDir, 'Lines', condNames{nCond}))
    end
    mkdir(fullfile(saveDir, 'Text', condNames{nCond}))
end

%% Repeated measures effects
%% MAIN EFFECT
% Lines
if numel(effectRM)>0
    if plotLines
        for nRm=1:numel(modalitiesRM)

            f=figure('units','centimeters','position',[0 0 10+4*numel(modalitiesRM{nRm}) 5+9/16*4*numel(modalitiesRM{nRm})],'visible','off');

            for x=1:numel(modalitiesRM{nRm})
                dataMeans(x)=nanmean(nanmean(data4plot.allData(:,col4means{nRm}(x,:)),2));
                dataSD(x)=nanstd(nanmean(data4plot.allData(:,col4means{nRm}(x,:)),2));
            end

            for x=1:numel(modalitiesRM{nRm})

                h=bar(x,dataMeans(x)); hold all
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

            nAov=findcolExact(aov.Effect,condNames{rmEffect(nRm)});
            pMAIN=aov{nAov,6};

            amp=[max(yl)-min(yl)];
            isSignificant=0;
            pValues=ones(1,size(postHoc.(condNames{rmEffect(nRm)}),1));
            pSelected=0;
            if pMAIN<pcritical(1) &  pSelected==0
                isSignificant=1;
                pValues=postHoc.(condNames{rmEffect(nRm)}){:,5};
                pSelected=1;
            end

            if pMAIN<pcritical(end) &  pSelected==0
                pValues=postHoc.(condNames{rmEffect(nRm)}){:,5};
                pSelected=1;
            end

            if pSelected==1
                nSignificant=1;
                for i=1:numel(pValues)
                    pV=pValues(i);
                    if pV<=pcritical(end)
                        if yl(2)>0
                            if pV<pcritical(1) & isSignificant==1
                                hline(yl(2)+0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNames{rmEffect(nRm)})(i,:),'lineWidth',1.5);
                            else
                                hline(yl(2)+0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNames{rmEffect(nRm)})(i,:),'lineWidth',1.5);
                            end
                        else
                            if pV<pcritical(1) & isSignificant==1
                                hline(yl(1)-0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNames{rmEffect(nRm)})(i,:),'lineWidth',1.5);
                            else
                                hline(yl(1)-0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNames{rmEffect(nRm)})(i,:),'lineWidth',1.5);
                            end
                        end
                        nSignificant=nSignificant+1;
                    end
                end
            end

            print('-dtiff',['-r' num2str(imageResolution)],fullfile(saveDir, 'Lines', condNames{rmEffect(nRm)}, 'All participants'))
            close
            clear ax yl dataMeans dataSD
        end
    end

    % Text
    for nRm=1:numel(effectRM)
        f=figure('units','centimeters','position',[0 0 10+4*numel(modalitiesRM{nRm}) 5+9/16*4*numel(modalitiesRM{nRm})],'visible','off');

        for x=1:numel(modalitiesRM{nRm})
            dataMeans(x)=nanmean(nanmean(data4plot.allData(:,col4means{nRm}(x,:)),2));
            dataSD(x)=nanstd(nanmean(data4plot.allData(:,col4means{nRm}(x,:)),2));
        end

        for x=1:numel(modalitiesRM{nRm})

            h=bar(x,dataMeans(x)); hold all
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

        nAov=findcolExact(aov.Effect,condNames{rmEffect(nRm)});
        pMAIN=aov{nAov,6};

        amp=[max(yl)-min(yl)];
        isSignificant=0;
        pValues=ones(1,size(postHoc.(condNames{rmEffect(nRm)}),1));
        pSelected=0;
        if pMAIN<pcritical(1) &  pSelected==0
            isSignificant=1;
            pValues=postHoc.(condNames{rmEffect(nRm)}){:,5};
            pSelected=1;
        end
        if pMAIN<pcritical(end) &  pSelected==0
            pValues=postHoc.(condNames{rmEffect(nRm)}){:,5};
            pSelected=1;
        end

        if plotSD==1
            dataMeans=dataMeans+sign(dataMeans).*dataSD;
        end

        nColSignificant=ones(1,numel(dataMeans));
        nSignificant=1;
        if pSelected==1
            for i=1:numel(pValues)
                pV=pValues(i);
                if pV<=pcritical(end)
                    if yl(2)>0
                        if abs(dataMeans(order4ES.(condNames{rmEffect(nRm)})(i,1)))>abs(dataMeans(order4ES.(condNames{rmEffect(nRm)})(i,2)))
                            addPvalue(order4ES.(condNames{rmEffect(nRm)})(i,1), dataMeans(order4ES.(condNames{rmEffect(nRm)})(i,1))+0.05*nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,1))*amp, pV, pcritical, colors{rmEffect(nRm)}(order4ES.(condNames{rmEffect(nRm)})(i,2),:))
                            nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,1))=nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,1))+1;
                        else
                            addPvalue(order4ES.(condNames{rmEffect(nRm)})(i,2), dataMeans(order4ES.(condNames{rmEffect(nRm)})(i,2))+0.05*nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,2))*amp, pV, pcritical, colors{rmEffect(nRm)}(order4ES.(condNames{rmEffect(nRm)})(i,1),:))
                            nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,2))=nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,2))+1;
                        end
                    else
                        if abs(dataMeans(order4ES.(condNames{rmEffect(nRm)})(i,1)))>abs(dataMeans(order4ES.(condNames{rmEffect(nRm)})(i,2)))
                            addPvalue(order4ES.(condNames{rmEffect(nRm)})(i,1), dataMeans(order4ES.(condNames{rmEffect(nRm)})(i,1))-0.05*nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,1))*amp, pV, pcritical, colors{rmEffect(nRm)}(order4ES.(condNames{rmEffect(nRm)})(i,2),:))
                            nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,1))=nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,1))+1;
                        else
                            addPvalue(order4ES.(condNames{rmEffect(nRm)})(i,2), dataMeans(order4ES.(condNames{rmEffect(nRm)})(i,2))-0.05*nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,2))*amp, pV, pcritical, colors{rmEffect(nRm)}(order4ES.(condNames{rmEffect(nRm)})(i,1),:))
                            nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,2))=nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,2))+1;
                        end
                    end
                    nSignificant=nSignificant+1;
                end
            end
        end

        yl=yl*1.05.^max(nColSignificant);
        set(ax,'ylim',[min(yl) max(yl)])

        print('-dtiff',['-r' num2str(imageResolution)],fullfile(saveDir, 'Text', condNames{rmEffect(nRm)}, 'All participants'))
        close
        clear ax yl dataMeans dataSD

    end
end

%% INTERACTIONS of RM
% lines
if numel(effectRM)>1
    if plotLines
        for nRm1=1:numel(modalitiesRM)
            for nRm2=1:numel(modalitiesRM)
                if nRm1~=nRm2

                    f=figure('units','centimeters','position',[0 0 10+4*numel(allModalities{rmEffect(nRm1)}) 5+numel(cond4effect{rmEffect(nRm2)})*numel(cond4effect{rmEffect(nRm1)})*9/16*4],'visible','off');

                    for nModRm=1:numel(cond4effect{rmEffect(nRm2)})

                        varData=tXl(:,col4means{nRm1}(:,nModRm));
                        dataMeans=nanmean(varData);
                        dataSD=nanstd(varData);

                        subplot(numel(cond4effect{rmEffect(nRm2)}),1,nModRm);
                        for x=1:numel(cond4effect{rmEffect(nRm1)})
                            h=bar(x,dataMeans(x)); hold all
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

                    nAovInt=findcolExact(aov.Effect,[condNames{rmEffect(nRm1)} ':' condNames{rmEffect(nRm2)}]);
                    nAov=findcolExact(aov.Effect,condNames{rmEffect(nRm1)});
                    pINT=aov{nAovInt,6};
                    pMAIN=aov{nAov,6};

                    amp=[max(max(yl))-min(min(yl))];
                    for nModRm=1:numel(cond4effect{rmEffect(nRm2)})

                        varData=tXl(:,col4means{nRm1}(:,nModRm));

                        isSignificant=0;
                        set(f,'CurrentAxes',ax{nModRm});
                        pValues=ones(1,size(postHoc.(condNames{rmEffect(nRm1)}),1));
                        pSelected=0;
                        if pINT<pcritical(1)
                            isSignificant=1;
                            phCut=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}]){:,1},condNames{rmEffect(nRm2)}{nModRm});
                            pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}]){phCut,6};
                            pSelected=1;
                        end
                        %                         if pMAIN<pcritical(1) &  pSelected==0
                        %                             isSignificant=1;
                        %                             pValues=postHoc.(condNames{rmEffect(nRm1)}){:,5};
                        %                             pSelected=1;
                        %                         end
                        if pINT<pcritical(end) &  pSelected==0
                            phCut=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}]){:,1},condNames{rmEffect(nRm2)}{nModRm});
                            pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}]){phCut,6};                        pSelected=1;
                        end
                        %                         if pMAIN<pcritical(end) &  pSelected==0
                        %                             pValues=postHoc.(condNames{rmEffect(nRm1)}){:,5};
                        %                             pSelected=1;
                        %                         end

                        if pSelected==1
                            nSignificant=1;
                            for i=1:numel(pValues)
                                pV=pValues(i);
                                if pV<=pcritical(end)
                                    if yl(2,nModRm)>0
                                        if pV<pcritical(1) & isSignificant==1
                                            hline(yl(2,nModRm)+0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNames{rmEffect(nRm1)})(i,:),'lineWidth',1.5);
                                        else
                                            hline(yl(2,nModRm)+0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNames{rmEffect(nRm1)})(i,:),'lineWidth',1.5);
                                        end
                                    else
                                        if pV<pcritical(1) & isSignificant==1
                                            hline(yl(1,nModRm)-0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNames{rmEffect(nRm1)})(i,:),'lineWidth',1.5);
                                        else
                                            hline(yl(1,nModRm)-0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNames{rmEffect(nRm1)})(i,:),'lineWidth',1.5);
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

                f=figure('units','centimeters','position',[0 0 10+4*numel(allModalities{rmEffect(nRm1)}) 5+numel(cond4effect{rmEffect(nRm2)})*numel(cond4effect{rmEffect(nRm1)})*9/16*4],'visible','off');
                for nModRm=1:numel(cond4effect{rmEffect(nRm2)})

                    varData=tXl(:,col4means{nRm1}(:,nModRm));
                    dataMeans=nanmean(varData);
                    dataSD=nanstd(varData);

                    subplot(numel(cond4effect{rmEffect(nRm2)}),1,nModRm);
                    for x=1:numel(cond4effect{rmEffect(nRm1)})
                        h=bar(x,dataMeans(x)); hold all
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

                nAovInt=findcolExact(aov.Effect,[condNames{rmEffect(nRm1)} ':' condNames{rmEffect(nRm2)}]);
                nAov=findcolExact(aov.Effect,condNames{rmEffect(nRm1)});
                pINT=aov{nAovInt,6};
                pMAIN=aov{nAov,6};

                amp=[max(max(yl))-min(min(yl))];
                for nModRm=1:numel(cond4effect{rmEffect(nRm2)})

                    varData=tXl(:,col4means{nRm1}(:,nModRm));
                    dataMeans=nanmean(varData);
                    dataSD=nanstd(varData);
                    if plotSD==1
                        dataMeans=dataMeans+sign(dataMeans).*dataSD;
                    end

                    isSignificant=0;
                    set(f,'CurrentAxes',ax{nModRm});
                    pValues=ones(1,size(postHoc.(condNames{rmEffect(nRm1)}),1));
                    pSelected=0;
                    if pINT<pcritical(1)
                        isSignificant=1;
                        phCut=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}]){:,1},condNames{rmEffect(nRm2)}{nModRm});
                        pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}]){phCut,6};
                        pSelected=1;
                    end
                    %                     if pMAIN<pcritical(1) &  pSelected==0
                    %                         isSignificant=1;
                    %                         pValues=postHoc.(condNames{rmEffect(nRm1)}){:,5};
                    %                         pSelected=1;
                    %                     end
                    if pINT<pcritical(end) &  pSelected==0
                        phCut=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}]){:,1},condNames{rmEffect(nRm2)}{nModRm});
                        pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}]){phCut,6};                        pSelected=1;
                    end
                    %                     if pMAIN<pcritical(end) &  pSelected==0
                    %                         pValues=postHoc.(condNames{rmEffect(nRm1)}){:,5};
                    %                         pSelected=1;
                    %                     end

                    nSignificant=1;
                    nColSignificant=ones(1,numel(dataMeans));
                    if pSelected==1
                        for i=1:numel(pValues)
                            pV=pValues(i);
                            if pV<=pcritical(end)
                                if yl(2)>0
                                    if abs(dataMeans(order4ES.(condNames{rmEffect(nRm1)})(i,1)))>abs(dataMeans(order4ES.(condNames{rmEffect(nRm1)})(i,2)))
                                        addPvalue(order4ES.(condNames{rmEffect(nRm1)})(i,1), dataMeans(order4ES.(condNames{rmEffect(nRm1)})(i,1))+0.022*numel(cond4effect{rmEffect(nRm2)})*nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,1))*amp, pV, pcritical, colors{rmEffect(nRm1)}(order4ES.(condNames{rmEffect(nRm1)})(i,2),:))
                                        nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,1))=nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,1))+1;

                                    else
                                        addPvalue(order4ES.(condNames{rmEffect(nRm1)})(i,2), dataMeans(order4ES.(condNames{rmEffect(nRm1)})(i,2))+0.022*numel(cond4effect{rmEffect(nRm2)})*nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,2))*amp, pV, pcritical, colors{rmEffect(nRm1)}(order4ES.(condNames{rmEffect(nRm1)})(i,1),:))
                                        nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,2))=nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,2))+1;
                                    end
                                else
                                    if abs(dataMeans(order4ES.(condNames{rmEffect(nRm1)})(i,1)))>abs(dataMeans(order4ES.(condNames{rmEffect(nRm1)})(i,2)))
                                        addPvalue(order4ES.(condNames{rmEffect(nRm1)})(i,1), dataMeans(order4ES.(condNames{rmEffect(nRm1)})(i,1))-0.022*numel(cond4effect{rmEffect(nRm2)})*nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,1))*amp, pV, pcritical, colors{rmEffect(nRm1)}(order4ES.(condNames{rmEffect(nRm1)})(i,2),:))
                                        nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,1))=nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,1))+1;
                                    else
                                        addPvalue(order4ES.(condNames{rmEffect(nRm1)})(i,2), dataMeans(order4ES.(condNames{rmEffect(nRm1)})(i,2))-0.022*numel(cond4effect{rmEffect(nRm2)})*nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,2))*amp, pV, pcritical, colors{rmEffect(nRm1)}(order4ES.(condNames{rmEffect(nRm1)})(i,1),:))
                                        nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,2))=nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,2))+1;
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

%% Main RM with simple IND interaction
% Lines
if numel(effectRM)>0 & numel(indEffect)>0
    if plotLines
        for nRm=1:numel(effectRM)
            for nInd=1:numel(indEffect)
                for nMod=1:numel(allMod.(condNames{indEffect(nInd)}))

                    f=figure('units','centimeters','position',[0 0 10+4*numel(modalitiesRM{nRm}) 5+9/16*4*numel(modalitiesRM{nRm})],'visible','off');

                    for x=1:numel(modalitiesRM{nRm})
                        dataMeans(x)=nanmean(nanmean(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){nMod})(:,col4means{nRm}(x,:)),2));
                        dataSD(x)=nanstd(nanmean(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){nMod})(:,col4means{nRm}(x,:)),2));
                    end

                    for x=1:numel(modalitiesRM{nRm})

                        h=bar(x,dataMeans(x)); hold all
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

                    nAovInt=findcolExact(aov.Effect,[condNames{indEffect(nInd)} ':' condNames{rmEffect(nRm)}]);
                    nAov=findcolExact(aov.Effect,condNames{rmEffect(nRm)});
                    pINT=aov{nAovInt,6};
                    pMAIN=aov{nAov,6};

                    amp=[max(yl)-min(yl)];
                    isSignificant=0;

                    pValues=ones(1,size(postHoc.(condNames{rmEffect(nRm)}),1));
                    pSelected=0;
                    if pINT<pcritical(1)
                        isSignificant=1;
                        rows=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd)}]){:,1},allMod.(condNames{indEffect(nInd)}){nMod});
                        pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd)}]){rows,6};
                        pSelected=1;
                    end
                    %                     if pMAIN<pcritical(1) &  pSelected==0
                    %                         isSignificant=1;
                    %                         pValues=postHoc.(condNames{rmEffect(nRm)}){:,5};
                    %                         pSelected=1;
                    %                     end
                    if pINT<pcritical(end) &  pSelected==0
                        rows=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd)}]){:,1},allMod.(condNames{indEffect(nInd)}){nMod});
                        pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd)}]){rows,6};
                        pSelected=1;
                    end
                    %                     if pMAIN<pcritical(end) &  pSelected==0
                    %                         pValues=postHoc.(condNames{rmEffect(nRm)}){:,5};
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
                                if yl(2)>0
                                    if pV<pcritical(1) & isSignificant==1
                                        hline(yl(2)+0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNames{rmEffect(nRm)})(i,:),'lineWidth',1.5);
                                    else
                                        hline(yl(2)+0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNames{rmEffect(nRm)})(i,:),'lineWidth',1.5);
                                    end
                                else
                                    if pV<pcritical(1) & isSignificant==1
                                        hline(yl(1)-0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNames{rmEffect(nRm)})(i,:),'lineWidth',1.5);
                                    else
                                        hline(yl(1)-0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNames{rmEffect(nRm)})(i,:),'lineWidth',1.5);
                                    end
                                end
                                nSignificant=nSignificant+1;
                            end
                        end
                    end


                    yl=yl*1.05.^max(nSignificant);
                    set(ax,'ylim',[min(yl) max(yl)])

                    print('-dtiff',['-r' num2str(imageResolution)], fullfile(saveDir, 'Lines', condNames{rmEffect(nRm)}, [condNames{indEffect(nInd)} ' = ' allMod.(condNames{indEffect(nInd)}){nMod}]))
                    close
                    clear ax yl dataMeans dataSD

                end
            end
        end
    end

    % Text
    for nRm=1:numel(effectRM)
        for nInd=1:numel(indEffect)
            for nMod=1:numel(allMod.(condNames{indEffect(nInd)}))

                f=figure('units','centimeters','position',[0 0 10+4*numel(modalitiesRM{nRm}) 5+9/16*4*numel(modalitiesRM{nRm})],'visible','off');

                for x=1:numel(modalitiesRM{nRm})
                    dataMeans(x)=nanmean(nanmean(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){nMod})(:,col4means{nRm}(x,:)),2));
                    dataSD(x)=nanstd(nanmean(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){nMod})(:,col4means{nRm}(x,:)),2));
                end

                for x=1:numel(modalitiesRM{nRm})

                    h=bar(x,dataMeans(x)); hold all
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

                nAovInt=findcolExact(aov.Effect,[condNames{indEffect(nInd)} ':' condNames{rmEffect(nRm)}]);
                nAov=findcolExact(aov.Effect,condNames{rmEffect(nRm)});
                pINT=aov{nAovInt,6};
                pMAIN=aov{nAov,6};

                amp=[max(yl)-min(yl)];
                isSignificant=0;

                pValues=ones(1,size(postHoc.(condNames{rmEffect(nRm)}),1));
                pSelected=0;
                if pINT<pcritical(1)
                    isSignificant=1;
                    rows=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd)}]){:,1},allMod.(condNames{indEffect(nInd)}){nMod});
                    pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd)}]){rows,6};
                    pSelected=1;
                end
                %                 if pMAIN<pcritical(1) &  pSelected==0
                %                     isSignificant=1;
                %                     pValues=postHoc.(condNames{rmEffect(nRm)}){:,5};
                %                     pSelected=1;
                %                 end
                if pINT<pcritical(end) &  pSelected==0
                    rows=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd)}]){:,1},allMod.(condNames{indEffect(nInd)}){nMod});
                    pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd)}]){rows,6};
                    pSelected=1;
                end
                %                 if pMAIN<pcritical(end) &  pSelected==0
                %                     pValues=postHoc.(condNames{rmEffect(nRm)}){:,5};
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
                            if yl(2)>0
                                if abs(dataMeans(order4ES.(condNames{rmEffect(nRm)})(i,1)))>abs(dataMeans(order4ES.(condNames{rmEffect(nRm)})(i,2)))
                                    addPvalue(order4ES.(condNames{rmEffect(nRm)})(i,1), dataMeans(order4ES.(condNames{rmEffect(nRm)})(i,1))+0.05*nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,1))*amp, pV, pcritical, colors{rmEffect(nRm)}(order4ES.(condNames{rmEffect(nRm)})(i,2),:))
                                    nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,1))=nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,1))+1;
                                else
                                    addPvalue(order4ES.(condNames{rmEffect(nRm)})(i,2), dataMeans(order4ES.(condNames{rmEffect(nRm)})(i,2))+0.05*nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,2))*amp, pV, pcritical, colors{rmEffect(nRm)}(order4ES.(condNames{rmEffect(nRm)})(i,1),:))
                                    nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,2))=nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,2))+1;

                                end
                            else
                                if abs(dataMeans(order4ES.(condNames{rmEffect(nRm)})(i,1)))>abs(dataMeans(order4ES.(condNames{rmEffect(nRm)})(i,2)))
                                    addPvalue(order4ES.(condNames{rmEffect(nRm)})(i,1), dataMeans(order4ES.(condNames{rmEffect(nRm)})(i,1))-0.05*nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,1))*amp, pV, pcritical, colors{rmEffect(nRm)}(order4ES.(condNames{rmEffect(nRm)})(i,2),:))
                                    nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,1))=nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,1))+1;
                                else
                                    addPvalue(order4ES.(condNames{rmEffect(nRm)})(i,2), dataMeans(order4ES.(condNames{rmEffect(nRm)})(i,2))-0.05*nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,2))*amp, pV, pcritical, colors{rmEffect(nRm)}(order4ES.(condNames{rmEffect(nRm)})(i,1),:))
                                    nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,2))=nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,2))+1;
                                end
                            end
                            nSignificant=nSignificant+1;
                        end
                    end
                end

                yl=yl*1.05.^max(nColSignificant);
                set(ax,'ylim',[min(yl) max(yl)])

                print('-dtiff',['-r' num2str(imageResolution)], fullfile(saveDir, 'Text', condNames{rmEffect(nRm)}, [condNames{indEffect(nInd)} ' = ' allMod.(condNames{indEffect(nInd)}){nMod}]))
                close
                clear ax yl dataMeans dataSD

            end
        end
    end
end

%% Main effect with double IND interaction
if numel(rmEffect)>0 & numel(indEffect)>1
    % Lines
    if plotLines
        for nRm=1:numel(effectRM)
            for nInd1=1:numel(indEffect)
                for nInd2=1:numel(indEffect)
                    if nInd2>nInd1

                        for nMod1=1:numel(allMod.(condNames{indEffect(nInd1)}))
                            for nMod2=1:numel(allMod.(condNames{indEffect(nInd2)}))

                                f=figure('units','centimeters','position',[0 0 10+4*numel(modalitiesRM{nRm}) 5+9/16*4*numel(modalitiesRM{nRm})],'visible','off');

                                for x=1:numel(modalitiesRM{nRm})
                                    dataMeans(x)=nanmean(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNames{indEffect(nInd1)}){nMod1}).(allMod.(condNames{indEffect(nInd2)}){nMod2})(:,col4means{nRm}(x,:)),2));
                                    dataSD(x)=nanstd(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNames{indEffect(nInd1)}){nMod1}).(allMod.(condNames{indEffect(nInd2)}){nMod2})(:,col4means{nRm}(x,:)),2));
                                end

                                for x=1:numel(modalitiesRM{nRm})

                                    h=bar(x,dataMeans(x)); hold all
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
                                nAovInt2=findcolExact(aov.Effect,[condNames{rmEffect(nRm)} ':' condNames{indEffect(nInd1)}]);
                                if isempty(nAovInt2)
                                    nAovInt2=findcolExact(aov.Effect,[condNames{rmEffect(nInd1)} ':' condNames{indEffect(nRm)}]);
                                end
                                nAovInt=findcolExact(aov.Effect,[condNames{rmEffect(nRm)} ':' condNames{indEffect(nInd2)}]);
                                if isempty(nAovInt)
                                    nAovInt=findcolExact(aov.Effect,[condNames{rmEffect(nInd2)} ':' condNames{indEffect(nRm)}]);
                                end
                                nAov=findcolExact(aov.Effect,condNames{rmEffect(nRm)});
                                pINT3=aov{nAovInt3,6};
                                pINT2=aov{nAovInt2,6};
                                pINT=aov{nAovInt,6};
                                pMAIN=aov{nAov,6};

                                amp=[max(yl)-min(yl)];
                                isSignificant=0;

                                pValues=ones(1,size(postHoc.(condNames{rmEffect(nRm)}),1));
                                pSelected=0;
                                if pINT3<pcritical(1)
                                    isSignificant=1;
                                    phCut1=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,1},allMod.(condNames{indEffect(nInd1)}){nMod1});
                                    phCut2=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,2},allMod.(condNames{indEffect(nInd2)}){nMod2});
                                    phCut=intersect(phCut1, phCut2);
                                    pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){phCut,7};
                                    pSelected=1;
                                end
                                %                                 if any([pINT<pcritical(1) pINT2<pcritical(1)]) &  pSelected==0
                                %                                     isSignificant=1;
                                %                                     pSelected=1;
                                %                                     if pINT<pINT2
                                %                                         phCut=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd2)}]){:,1}, allMod.(condNames{indEffect(nInd2)}){nMod2});
                                %                                         pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd2)}]){phCut,6};
                                %                                     else
                                %                                         phCut=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)}]){:,1}, allMod.(condNames{indEffect(nInd1)}){nMod1});
                                %                                         pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)}]){phCut,6};
                                %                                     end
                                %                                 end
                                %                                 if pMAIN<pcritical(1) &  pSelected==0
                                %                                     isSignificant=1;
                                %                                     pValues=postHoc.(condNames{rmEffect(nRm)}){:,5};
                                %                                     pSelected=1;
                                %                                 end
                                if pINT3<pcritical(end) &  pSelected==0
                                    phCut1=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,1},allMod.(condNames{indEffect(nInd1)}){nMod1});
                                    phCut2=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,2},allMod.(condNames{indEffect(nInd2)}){nMod2});
                                    phCut=intersect(phCut1, phCut2);
                                    pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){phCut,7};
                                    pSelected=1;
                                end
                                %                                 if any([pINT<pcritical(end) pINT2<pcritical(end)]) &  pSelected==0
                                %                                     pSelected=1;
                                %                                     if pINT<pINT2
                                %                                         phCut=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd2)}]){:,1}, allMod.(condNames{indEffect(nInd2)}){nMod2});
                                %                                         pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd2)}]){phCut,6};
                                %                                     else
                                %                                         phCut=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)}]){:,1}, allMod.(condNames{indEffect(nInd1)}){nMod1});
                                %                                         pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)}]){phCut,6};
                                %                                     end
                                %                                 end
                                %                                 if pMAIN<pcritical(end) &  pSelected==0
                                %                                     pValues=postHoc.(condNames{rmEffect(nRm)}){:,5};
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
                                            if yl(2)>0
                                                if pV<pcritical(1) & isSignificant==1
                                                    hline(yl(2)+0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNames{rmEffect(nRm)})(i,:),'lineWidth',1.5);
                                                else
                                                    hline(yl(2)+0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNames{rmEffect(nRm)})(i,:),'lineWidth',1.5);
                                                end
                                            else
                                                if pV<pcritical(1) & isSignificant==1
                                                    hline(yl(1)-0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNames{rmEffect(nRm)})(i,:),'lineWidth',1.5);
                                                else
                                                    hline(yl(1)-0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNames{rmEffect(nRm)})(i,:),'lineWidth',1.5);
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

                    for nMod1=1:numel(allMod.(condNames{indEffect(nInd1)}))
                        for nMod2=1:numel(allMod.(condNames{indEffect(nInd2)}))

                            f=figure('units','centimeters','position',[0 0 10+4*numel(modalitiesRM{nRm}) 5+9/16*4*numel(modalitiesRM{nRm})],'visible','off');

                            for x=1:numel(modalitiesRM{nRm})
                                dataMeans(x)=nanmean(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNames{indEffect(nInd1)}){nMod1}).(allMod.(condNames{indEffect(nInd2)}){nMod2})(:,col4means{nRm}(x,:)),2));
                                dataSD(x)=nanstd(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNames{indEffect(nInd1)}){nMod1}).(allMod.(condNames{indEffect(nInd2)}){nMod2})(:,col4means{nRm}(x,:)),2));
                            end

                            for x=1:numel(modalitiesRM{nRm})

                                h=bar(x,dataMeans(x)); hold all
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
                            %                             nAov=findcolExact(aov.Effect,condNames{rmEffect(nRm)});
                            pINT3=aov{nAovInt3,6};
                            %                             pINT2=aov{nAovInt2,6};
                            %                             pINT=aov{nAovInt,6};
                            %                             pMAIN=aov{nAov,6};

                            amp=[max(yl)-min(yl)];
                            isSignificant=0;

                            pValues=ones(1,size(postHoc.(condNames{rmEffect(nRm)}),1));
                            pSelected=0;
                            if pINT3<pcritical(1)
                                isSignificant=1;
                                phCut1=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,1},allMod.(condNames{indEffect(nInd1)}){nMod1});
                                phCut2=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,2},allMod.(condNames{indEffect(nInd2)}){nMod2});
                                phCut=intersect(phCut1, phCut2);
                                pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){phCut,7};
                                pSelected=1;
                            end
                            %                             if any([pINT<pcritical(1) pINT2<pcritical(1)]) &  pSelected==0
                            %                                 isSignificant=1;
                            %                                 pSelected=1;
                            %                                 if pINT<pINT2
                            %                                     phCut=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd2)}]){:,1}, allMod.(condNames{indEffect(nInd2)}){nMod2});
                            %                                     pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd2)}]){phCut,6};
                            %                                 else
                            %                                     phCut=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)}]){:,1}, allMod.(condNames{indEffect(nInd1)}){nMod1});
                            %                                     pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)}]){phCut,6};
                            %                                 end
                            %                             end
                            %                             if pMAIN<pcritical(1) &  pSelected==0
                            %                                 isSignificant=1;
                            %                                 pValues=postHoc.(condNames{rmEffect(nRm)}){:,5};
                            %                                 pSelected=1;
                            %                             end
                            if pINT3<pcritical(end) &  pSelected==0
                                phCut1=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,1},allMod.(condNames{indEffect(nInd1)}){nMod1});
                                phCut2=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,2},allMod.(condNames{indEffect(nInd2)}){nMod2});
                                phCut=intersect(phCut1, phCut2);
                                pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){phCut,7};
                                pSelected=1;
                            end
                            %                             if any([pINT<pcritical(end) pINT2<pcritical(end)]) &  pSelected==0
                            %                                 pSelected=1;
                            %                                 if pINT<pINT2
                            %                                     phCut=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd2)}]){:,1}, allMod.(condNames{indEffect(nInd2)}){nMod2});
                            %                                     pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd2)}]){phCut,6};
                            %                                 else
                            %                                     phCut=findcolExact(postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)}]){:,1}, allMod.(condNames{indEffect(nInd1)}){nMod1});
                            %                                     pValues=postHoc.([condNames{rmEffect(nRm)} 'By' condNames{indEffect(nInd1)}]){phCut,6};
                            %                                 end
                            %                             end
                            %                             if pMAIN<pcritical(end) &  pSelected==0
                            %                                 pValues=postHoc.(condNames{rmEffect(nRm)}){:,5};
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
                                        if yl(2)>0
                                            if abs(dataMeans(order4ES.(condNames{rmEffect(nRm)})(i,1)))>abs(dataMeans(order4ES.(condNames{rmEffect(nRm)})(i,2)))
                                                addPvalue(order4ES.(condNames{rmEffect(nRm)})(i,1), dataMeans(order4ES.(condNames{rmEffect(nRm)})(i,1))+0.05*nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,1))*amp, pV, pcritical, colors{rmEffect(nRm)}(order4ES.(condNames{rmEffect(nRm)})(i,2),:))
                                                nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,1))=nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,1))+1;
                                            else
                                                addPvalue(order4ES.(condNames{rmEffect(nRm)})(i,2), dataMeans(order4ES.(condNames{rmEffect(nRm)})(i,2))+0.05*nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,2))*amp, pV, pcritical, colors{rmEffect(nRm)}(order4ES.(condNames{rmEffect(nRm)})(i,1),:))
                                                nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,2))=nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,2))+1;

                                            end
                                        else
                                            if abs(dataMeans(order4ES.(condNames{rmEffect(nRm)})(i,1)))>abs(dataMeans(order4ES.(condNames{rmEffect(nRm)})(i,2)))
                                                addPvalue(order4ES.(condNames{rmEffect(nRm)})(i,1), dataMeans(order4ES.(condNames{rmEffect(nRm)})(i,1))-0.05*nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,1))*amp, pV, pcritical, colors{rmEffect(nRm)}(order4ES.(condNames{rmEffect(nRm)})(i,2),:))
                                                nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,1))=nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,1))+1;
                                            else
                                                addPvalue(order4ES.(condNames{rmEffect(nRm)})(i,2), dataMeans(order4ES.(condNames{rmEffect(nRm)})(i,2))-0.05*nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,2))*amp, pV, pcritical, colors{rmEffect(nRm)}(order4ES.(condNames{rmEffect(nRm)})(i,1),:))
                                                nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,2))=nColSignificant(order4ES.(condNames{rmEffect(nRm)})(i,2))+1;
                                            end
                                        end
                                        nSignificant=nSignificant+1;
                                    end
                                end
                            end

                            yl=yl*1.05.^max(nColSignificant);
                            set(ax,'ylim',[min(yl) max(yl)])

                            print('-dtiff',['-r' num2str(imageResolution)], fullfile(saveDir, 'Text', condNames{rmEffect(nRm)}, [condNames{indEffect(nInd1)} ' By ' condNames{indEffect(nInd2)} ' = ' allMod.(condNames{indEffect(nInd1)}){nMod1} ' ' allMod.(condNames{indEffect(nInd2)}){nMod2}]))
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
    if plotLines
        for nRm1=1:numel(modalitiesRM)
            for nRm2=1:numel(modalitiesRM)
                if nRm1~=nRm2
                    for nInd=1:numel(indEffect)
                        for nMod=1:numel(allMod.(condNames{indEffect(nInd)}))

                            f=figure('units','centimeters','position',[0 0 10+4*numel(allModalities{rmEffect(nRm1)}) 5+numel(cond4effect{rmEffect(nRm2)})*numel(cond4effect{rmEffect(nRm1)})*9/16*4],'visible','off');

                            for nModRm=1:numel(cond4effect{rmEffect(nRm2)})
                                dataMeans=nanmean(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){nMod})(:,col4means{nRm1}(:,nModRm)));
                                dataSD=nanstd(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){nMod})(:,col4means{nRm1}(:,nModRm)));

                                subplot(numel(cond4effect{rmEffect(nRm2)}),1,nModRm);
                                for x=1:numel(cond4effect{rmEffect(nRm1)})
                                    h=bar(x,dataMeans(x)); hold all
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

                            nAov=findcolExact(aov.Effect,condNames{rmEffect(nRm1)});
                            pINT3=aov{nAovInt3,6};
                            pINT2=aov{nAovInt2,6};
                            pINT=aov{nAovInt,6};
                            pMAIN=aov{nAov,6};

                            amp=[max(max(yl))-min(min(yl))];
                            for nModRm=1:numel(cond4effect{rmEffect(nRm2)})
                                isSignificant=0;
                                set(f,'CurrentAxes',ax{nModRm});
                                pValues=ones(1,size(postHoc.(condNames{rmEffect(nRm1)}),1));
                                pSelected=0;

                                if pINT3<pcritical(1)
                                    isSignificant=1;
                                    phCut1=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd)}]){:,1}, allMod.(condNames{indEffect(nInd)}){nMod});
                                    phCut2=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd)}]){:,2}, allMod.(condNames{rmEffect(nRm2)}){nModRm});
                                    phCut=intersect(phCut1, phCut2);
                                    pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd)}]){phCut,7};
                                    pSelected=1;
                                end
                                %                                 if any([pINT<pcritical(1) pINT2<pcritical(1)]) &  pSelected==0
                                %                                     isSignificant=1;
                                %                                     pSelected=1;
                                %                                     if pINT<pINT2
                                %                                         phCut=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}]){:,1}, allMod.(condNames{rmEffect(nRm2)}){nModRm});
                                %                                         pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}]){phCut,6};
                                %                                     else
                                %                                         phCut=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{indEffect(nInd)}]){:,1}, allMod.(condNames{indEffect(nInd)}){nMod});
                                %                                         pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{indEffect(nInd)}]){phCut,6};
                                %                                     end
                                %                                 end
                                %                                 if pMAIN<pcritical(1) &  pSelected==0
                                %                                     isSignificant=1;
                                %                                     pValues=postHoc.(condNames{rmEffect(nRm1)}){:,5};
                                %                                     pSelected=1;
                                %                                 end
                                if pINT3<pcritical(end) &  pSelected==0
                                    phCut1=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd)}]){:,1}, allMod.(condNames{indEffect(nInd)}){nMod});
                                    phCut2=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd)}]){:,2}, allMod.(condNames{rmEffect(nRm2)}){nModRm});
                                    phCut=intersect(phCut1, phCut2);
                                    pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd)}]){phCut,7};
                                    pSelected=1;
                                end
                                %                                 if any([pINT<pcritical(end) pINT2<pcritical(end)]) &  pSelected==0
                                %                                     pSelected=1;
                                %                                     if pINT<pINT2
                                %                                         phCut=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}]){:,1}, allMod.(condNames{rmEffect(nRm2)}){nModRm});
                                %                                         pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}]){phCut,6};
                                %                                     else
                                %                                         phCut=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{indEffect(nInd)}]){:,1}, allMod.(condNames{indEffect(nInd)}){nMod});
                                %                                         pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{indEffect(nInd)}]){phCut,6};
                                %                                     end
                                %                                 end
                                %                                 if pMAIN<pcritical(end) &  pSelected==0
                                %                                     pValues=postHoc.(condNames{rmEffect(nRm1)}){:,5};
                                %                                     pSelected=1;
                                %                                 end

                                if pSelected==1
                                    nSignificant=1;
                                    for i=1:numel(pValues)
                                        pV=pValues(i);
                                        if pV<=pcritical(end)
                                            if yl(2,nModRm)>0
                                                if pV<pcritical(1) & isSignificant==1
                                                    hline(yl(2,nModRm)+0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNames{rmEffect(nRm1)})(i,:),'lineWidth',1.5);
                                                else
                                                    hline(yl(2,nModRm)+0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNames{rmEffect(nRm1)})(i,:),'lineWidth',1.5);
                                                end
                                            else
                                                if pV<pcritical(1) & isSignificant==1
                                                    hline(yl(1,nModRm)-0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNames{rmEffect(nRm1)})(i,:),'lineWidth',1.5);
                                                else
                                                    hline(yl(1,nModRm)-0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNames{rmEffect(nRm1)})(i,:),'lineWidth',1.5);
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

                            print('-dtiff',['-r' num2str(imageResolution)],fullfile(saveDir, 'Lines', condNames{rmEffect(nRm1)}, [condNames{indEffect(nInd)} ' = ' allMod.(condNames{indEffect(nInd)}){nMod} ' By ' condNames{rmEffect(nRm2)}]))
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
                    for nMod=1:numel(allMod.(condNames{indEffect(nInd)}))

                        f=figure('units','centimeters','position',[0 0 10+4*numel(allModalities{rmEffect(nRm1)}) 5+numel(cond4effect{rmEffect(nRm2)})*numel(cond4effect{rmEffect(nRm1)})*9/16*4],'visible','off');

                        for nModRm=1:numel(cond4effect{rmEffect(nRm2)})
                            dataMeans=nanmean(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){nMod})(:,col4means{nRm1}(:,nModRm)));
                            dataSD=nanstd(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){nMod})(:,col4means{nRm1}(:,nModRm)));

                            subplot(numel(cond4effect{rmEffect(nRm2)}),1,nModRm);

                            for x=1:numel(cond4effect{rmEffect(nRm1)})
                                h=bar(x,dataMeans(x)); hold all
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
                        nAov=findcolExact(aov.Effect,condNames{rmEffect(nRm1)});
                        pINT3=aov{nAovInt3,6};
                        pINT2=aov{nAovInt2,6};
                        pINT=aov{nAovInt,6};
                        pMAIN=aov{nAov,6};

                        amp=[max(max(yl))-min(min(yl))];
                        for nModRm=1:numel(cond4effect{rmEffect(nRm2)})

                            dataMeans=nanmean(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){nMod})(:,col4means{nRm1}(:,nModRm)));
                            dataSD=nanstd(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){nMod})(:,col4means{nRm1}(:,nModRm)));

                            if plotSD==1
                                dataMeans=dataMeans+sign(dataMeans).*dataSD;
                            end

                            isSignificant=0;
                            set(f,'CurrentAxes',ax{nModRm});
                            pValues=ones(1,size(postHoc.(condNames{rmEffect(nRm1)}),1));
                            pSelected=0;
                            if pINT3<pcritical(1)
                                isSignificant=1;
                                phCut1=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd)}]){:,1}, allMod.(condNames{indEffect(nInd)}){nMod});
                                phCut2=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd)}]){:,2}, allMod.(condNames{rmEffect(nRm2)}){nModRm});
                                phCut=intersect(phCut1, phCut2);
                                pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd)}]){phCut,7};
                                pSelected=1;
                            end
                            %                             if any([pINT<pcritical(1) pINT2<pcritical(1)]) &  pSelected==0
                            %                                 isSignificant=1;
                            %                                 pSelected=1;
                            %                                 if pINT<pINT2
                            %                                     phCut=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}]){:,1}, allMod.(condNames{rmEffect(nRm2)}){nModRm});
                            %                                     pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}]){phCut,6};
                            %                                 else
                            %                                     phCut=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{indEffect(nInd)}]){:,1}, allMod.(condNames{indEffect(nInd)}){nMod});
                            %                                     pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{indEffect(nInd)}]){phCut,6};
                            %                                 end
                            %                             end
                            %                             if pMAIN<pcritical(1) &  pSelected==0
                            %                                 isSignificant=1;
                            %                                 pValues=postHoc.(condNames{rmEffect(nRm1)}){:,5};
                            %                                 pSelected=1;
                            %                             end
                            if pINT3<pcritical(end) &  pSelected==0
                                phCut1=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd)}]){:,1}, allMod.(condNames{indEffect(nInd)}){nMod});
                                phCut2=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd)}]){:,2}, allMod.(condNames{rmEffect(nRm2)}){nModRm});
                                phCut=intersect(phCut1, phCut2);
                                pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd)}]){phCut,7};
                                pSelected=1;
                            end
                            %                             if any([pINT<pcritical(end) pINT2<pcritical(end)]) &  pSelected==0
                            %                                 pSelected=1;
                            %                                 if pINT<pINT2
                            %                                     phCut=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}]){:,1}, allMod.(condNames{rmEffect(nRm2)}){nModRm});
                            %                                     pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)}]){phCut,6};
                            %                                 else
                            %                                     phCut=findcol(postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{indEffect(nInd)}]){:,1}, allMod.(condNames{indEffect(nInd)}){nMod});
                            %                                     pValues=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{indEffect(nInd)}]){phCut,6};
                            %                                 end
                            %                             end
                            %                             if pMAIN<pcritical(end) &  pSelected==0
                            %                                 pValues=postHoc.(condNames{rmEffect(nRm1)}){:,5};
                            %                                 pSelected=1;
                            %                             end

                            nSignificant=1;
                            nColSignificant=ones(1,numel(dataMeans));
                            if pSelected==1
                                for i=1:numel(pValues)
                                    pV=pValues(i);
                                    if pV<=pcritical(end)
                                        if yl(2)>0
                                            if abs(dataMeans(order4ES.(condNames{rmEffect(nRm1)})(i,1)))>abs(dataMeans(order4ES.(condNames{rmEffect(nRm1)})(i,2)))
                                                addPvalue(order4ES.(condNames{rmEffect(nRm1)})(i,1), dataMeans(order4ES.(condNames{rmEffect(nRm1)})(i,1))+0.022*numel(cond4effect{rmEffect(nRm2)})*nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,1))*amp, pV, pcritical, colors{rmEffect(nRm1)}(order4ES.(condNames{rmEffect(nRm1)})(i,2),:))
                                                nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,1))=nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,1))+1;

                                            else
                                                addPvalue(order4ES.(condNames{rmEffect(nRm1)})(i,2), dataMeans(order4ES.(condNames{rmEffect(nRm1)})(i,2))+0.022*numel(cond4effect{rmEffect(nRm2)})*nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,2))*amp, pV, pcritical, colors{rmEffect(nRm1)}(order4ES.(condNames{rmEffect(nRm1)})(i,1),:))
                                                nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,2))=nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,2))+1;
                                            end
                                        else
                                            if abs(dataMeans(order4ES.(condNames{rmEffect(nRm1)})(i,1)))>abs(dataMeans(order4ES.(condNames{rmEffect(nRm1)})(i,2)))
                                                addPvalue(order4ES.(condNames{rmEffect(nRm1)})(i,1), dataMeans(order4ES.(condNames{rmEffect(nRm1)})(i,1))-0.022*numel(cond4effect{rmEffect(nRm2)})*nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,1))*amp, pV, pcritical, colors{rmEffect(nRm1)}(order4ES.(condNames{rmEffect(nRm1)})(i,2),:))
                                                nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,1))=nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,1))+1;
                                            else
                                                addPvalue(order4ES.(condNames{rmEffect(nRm1)})(i,2), dataMeans(order4ES.(condNames{rmEffect(nRm1)})(i,2))-0.022*numel(cond4effect{rmEffect(nRm2)})*nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,2))*amp, pV, pcritical, colors{rmEffect(nRm1)}(order4ES.(condNames{rmEffect(nRm1)})(i,1),:))
                                                nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,2))=nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,2))+1;
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

                        print('-dtiff',['-r' num2str(imageResolution)],fullfile(saveDir, 'Text', condNames{rmEffect(nRm1)}, [condNames{indEffect(nInd)} ' = ' allMod.(condNames{indEffect(nInd)}){nMod} ' By ' condNames{rmEffect(nRm2)}]))
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
    if plotLines
        for nRm1=1:numel(modalitiesRM)
            for nRm2=1:numel(modalitiesRM)
                if nRm1~=nRm2
                    for nInd1=1:numel(indEffect)
                        for nInd2=1:numel(indEffect)
                            if nInd2>nInd1
                                for nModInd1=1:numel(allMod.(condNames{indEffect(nInd1)}))
                                    for nModInd2=1:numel(allMod.(condNames{indEffect(nInd2)}))

                                        f=figure('units','centimeters','position',[0 0 10+4*numel(allModalities{rmEffect(nRm1)}) 5+numel(cond4effect{rmEffect(nRm2)})*numel(cond4effect{rmEffect(nRm1)})*9/16*4],'visible','off');

                                        for nModRm=1:numel(cond4effect{rmEffect(nRm2)})
                                            dataMeans=nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNames{indEffect(nInd1)}){nModInd1}).(allMod.(condNames{indEffect(nInd2)}){nModInd2})(:,col4means{nRm1}(:,nModRm)));
                                            dataSD=nanstd(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNames{indEffect(nInd1)}){nModInd1}).(allMod.(condNames{indEffect(nInd2)}){nModInd2})(:,col4means{nRm1}(:,nModRm)));

                                            subplot(numel(cond4effect{rmEffect(nRm2)}),1,nModRm);
                                            for x=1:numel(cond4effect{rmEffect(nRm1)})
                                                h=bar(x,dataMeans(x)); hold all
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

                                            dataMeans=nanmean(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){nMod})(:,col4means{nRm1}(:,nModRm)));
                                            dataSD=nanstd(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){nMod})(:,col4means{nRm1}(:,nModRm)));

                                            if plotSD==1
                                                dataMeans=dataMeans+sign(dataMeans).*dataSD;
                                            end

                                            isSignificant=0;
                                            set(f,'CurrentAxes',ax{nModRm});
                                            pValues=ones(1,size(postHoc.(condNames{rmEffect(nRm1)}),1));
                                            pSelected=0;
                                            if pINT4<pcritical(1)
                                                isSignificant=1;
                                                ph4=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]);
                                                phCut1=findcol(ph4{:,1}, allMod.(condNames{indEffect(nInd1)}){nModInd1});
                                                phCut2=findcol(ph4{:,2}, allMod.(condNames{indEffect(nInd2)}){nModInd2});
                                                phCut3=findcol(ph4{:,3}, allMod.(condNames{rmEffect(nRm2)}){nModRm});
                                                phCut4=intersect(phCut1, phCut2);
                                                phCut=intersect(phCut4, phCut3);
                                                pValues=ph4{phCut,8};
                                                pSelected=1;
                                            end
                                            if pINT4<pcritical(end) & pSelected==0
                                                ph4=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]);
                                                phCut1=findcol(ph4{:,1}, allMod.(condNames{indEffect(nInd1)}){nModInd1});
                                                phCut2=findcol(ph4{:,2}, allMod.(condNames{indEffect(nInd2)}){nModInd2});
                                                phCut3=findcol(ph4{:,3}, allMod.(condNames{rmEffect(nRm2)}){nModRm});
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
                                                                hline(yl(2,nModRm)+0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNames{rmEffect(nRm1)})(i,:),'lineWidth',1.5);
                                                            else
                                                                hline(yl(2,nModRm)+0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNames{rmEffect(nRm1)})(i,:),'lineWidth',1.5);
                                                            end
                                                        else
                                                            if pV<pcritical(1) & isSignificant==1
                                                                hline(yl(1,nModRm)-0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNames{rmEffect(nRm1)})(i,:),'lineWidth',1.5);
                                                            else
                                                                hline(yl(1,nModRm)-0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNames{rmEffect(nRm1)})(i,:),'lineWidth',1.5);
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

                                        print('-dtiff',['-r' num2str(imageResolution)],fullfile(saveDir, 'Lines', condNames{rmEffect(nRm1)}, [condNames{indEffect(nInd1)} ' By ' condNames{indEffect(nInd2)} ' = ' allMod.(condNames{indEffect(nInd1)}){nModInd1} ' ' allMod.(condNames{indEffect(nInd2)}){nModInd2} ' By ' condNames{rmEffect(nRm2)}]))
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
                            for nModInd1=1:numel(allMod.(condNames{indEffect(nInd1)}))
                                for nModInd2=1:numel(allMod.(condNames{indEffect(nInd2)}))

                                    f=figure('units','centimeters','position',[0 0 10+4*numel(allModalities{rmEffect(nRm1)}) 5+numel(cond4effect{rmEffect(nRm2)})*numel(cond4effect{rmEffect(nRm1)})*9/16*4],'visible','off');

                                    for nModRm=1:numel(cond4effect{rmEffect(nRm2)})
                                        dataMeans=nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNames{indEffect(nInd1)}){nModInd1}).(allMod.(condNames{indEffect(nInd2)}){nModInd2})(:,col4means{nRm1}(:,nModRm)));
                                        dataSD=nanstd(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNames{indEffect(nInd1)}){nModInd1}).(allMod.(condNames{indEffect(nInd2)}){nModInd2})(:,col4means{nRm1}(:,nModRm)));

                                        subplot(numel(cond4effect{rmEffect(nRm2)}),1,nModRm);
                                        for x=1:numel(cond4effect{rmEffect(nRm1)})
                                            h=bar(x,dataMeans(x)); hold all
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

                                        dataMeans=nanmean(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){nMod})(:,col4means{nRm1}(:,nModRm)));
                                        dataSD=nanstd(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){nMod})(:,col4means{nRm1}(:,nModRm)));

                                        if plotSD==1
                                            dataMeans=dataMeans+sign(dataMeans).*dataSD;
                                        end

                                        isSignificant=0;
                                        set(f,'CurrentAxes',ax{nModRm});
                                        pValues=ones(1,size(postHoc.(condNames{rmEffect(nRm1)}),1));
                                        pSelected=0;
                                        if pINT4<pcritical(1)
                                            isSignificant=1;
                                            ph4=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]);
                                            phCut1=findcol(ph4{:,1}, allMod.(condNames{indEffect(nInd1)}){nModInd1});
                                            phCut2=findcol(ph4{:,2}, allMod.(condNames{indEffect(nInd2)}){nModInd2});
                                            phCut3=findcol(ph4{:,3}, allMod.(condNames{rmEffect(nRm2)}){nModRm});
                                            phCut4=intersect(phCut1, phCut2);
                                            phCut=intersect(phCut4, phCut3);
                                            pValues=ph4{phCut,8};
                                            pSelected=1;
                                        end
                                        if pINT4<pcritical(end) & pSelected==0
                                            ph4=postHoc.([condNames{rmEffect(nRm1)} 'By' condNames{rmEffect(nRm2)} 'By' condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]);
                                            phCut1=findcol(ph4{:,1}, allMod.(condNames{indEffect(nInd1)}){nModInd1});
                                            phCut2=findcol(ph4{:,2}, allMod.(condNames{indEffect(nInd2)}){nModInd2});
                                            phCut3=findcol(ph4{:,3}, allMod.(condNames{rmEffect(nRm2)}){nModRm});
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
                                                    if yl(2)>0
                                                        if abs(dataMeans(order4ES.(condNames{rmEffect(nRm1)})(i,1)))>abs(dataMeans(order4ES.(condNames{rmEffect(nRm1)})(i,2)))
                                                            addPvalue(order4ES.(condNames{rmEffect(nRm1)})(i,1), dataMeans(order4ES.(condNames{rmEffect(nRm1)})(i,1))+0.022*numel(cond4effect{rmEffect(nRm2)})*nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,1))*amp, pV, pcritical, colors{rmEffect(nRm1)}(order4ES.(condNames{rmEffect(nRm1)})(i,2),:))
                                                            nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,1))=nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,1))+1;

                                                        else
                                                            addPvalue(order4ES.(condNames{rmEffect(nRm1)})(i,2), dataMeans(order4ES.(condNames{rmEffect(nRm1)})(i,2))+0.022*numel(cond4effect{rmEffect(nRm2)})*nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,2))*amp, pV, pcritical, colors{rmEffect(nRm1)}(order4ES.(condNames{rmEffect(nRm1)})(i,1),:))
                                                            nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,2))=nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,2))+1;
                                                        end
                                                    else
                                                        if abs(dataMeans(order4ES.(condNames{rmEffect(nRm1)})(i,1)))>abs(dataMeans(order4ES.(condNames{rmEffect(nRm1)})(i,2)))
                                                            addPvalue(order4ES.(condNames{rmEffect(nRm1)})(i,1), dataMeans(order4ES.(condNames{rmEffect(nRm1)})(i,1))-0.022*numel(cond4effect{rmEffect(nRm2)})*nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,1))*amp, pV, pcritical, colors{rmEffect(nRm1)}(order4ES.(condNames{rmEffect(nRm1)})(i,2),:))
                                                            nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,1))=nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,1))+1;
                                                        else
                                                            addPvalue(order4ES.(condNames{rmEffect(nRm1)})(i,2), dataMeans(order4ES.(condNames{rmEffect(nRm1)})(i,2))-0.022*numel(cond4effect{rmEffect(nRm2)})*nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,2))*amp, pV, pcritical, colors{rmEffect(nRm1)}(order4ES.(condNames{rmEffect(nRm1)})(i,1),:))
                                                            nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,2))=nColSignificant(order4ES.(condNames{rmEffect(nRm1)})(i,2))+1;
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

                                    print('-dtiff',['-r' num2str(imageResolution)],fullfile(saveDir, 'Text', condNames{rmEffect(nRm1)}, [condNames{indEffect(nInd1)} ' By ' condNames{indEffect(nInd2)} ' = ' allMod.(condNames{indEffect(nInd1)}){nModInd1} ' ' allMod.(condNames{indEffect(nInd2)}){nModInd2} ' By ' condNames{rmEffect(nRm2)}]))
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
    if plotLines
        for nInd=1:numel(indEffect)

            f=figure('units','centimeters','position',[0 0 10+4*numel(modalitiesInd{nInd}) 5+9/16*4*numel(modalitiesInd{nInd})],'visible','off');

            for nMod=1:numel(modalitiesInd{nInd})
                dataMeans(nMod)=nanmean(nanmean(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){nMod}),2));
                dataSD(nMod)=nanstd(nanmean(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){nMod}),2));
            end

            for x=1:numel(modalitiesInd{nInd})

                h=bar(x,dataMeans(x)); hold all
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
            xlabel(condNames{indEffect(nInd)})
            ylabel(units)
            yl=ylim;
            box off
            ax=gca;
            ax.XGrid='off';
            ax.YGrid='on';

            nAov=findcolExact(aov.Effect,condNames{indEffect(nInd)});
            pMAIN=aov{nAov,6};

            amp=[max(yl)-min(yl)];
            isSignificant=0;
            pValues=ones(1,size(postHoc.(condNames{indEffect(nInd)}),1));
            pSelected=0;
            if pMAIN<pcritical(1) &  pSelected==0
                isSignificant=1;
                pValues=postHoc.(condNames{indEffect(nInd)}){:,5};
                pSelected=1;
            end

            if pMAIN<pcritical(end) &  pSelected==0
                pValues=postHoc.(condNames{indEffect(nInd)}){:,5};
                pSelected=1;
            end

            nSignificant=1;
            if pSelected==1
                nSignificant=1;
                for i=1:numel(pValues)
                    pV=pValues(i);
                    if pV<=pcritical(end)
                        if yl(2)>0
                            if pV<pcritical(1) & isSignificant==1
                                hline(yl(2)+0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNames{indEffect(nInd)})(i,:),'lineWidth',1.5);
                            else
                                hline(yl(2)+0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNames{indEffect(nInd)})(i,:),'lineWidth',1.5);
                            end
                        else
                            if pV<pcritical(1) & isSignificant==1
                                hline(yl(1)-0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNames{indEffect(nInd)})(i,:),'lineWidth',1.5);
                            else
                                hline(yl(1)-0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNames{indEffect(nInd)})(i,:),'lineWidth',1.5);
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
            clear ax yl dataMeans dataSD

        end
    end

    % Text
    for nInd=1:numel(indEffect)

        f=figure('units','centimeters','position',[0 0 10+4*numel(modalitiesInd{nInd}) 5+9/16*4*numel(modalitiesInd{nInd})],'visible','off');

        for nMod=1:numel(modalitiesInd{nInd})
            dataMeans(nMod)=nanmean(nanmean(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){nMod}),2));
            dataSD(nMod)=nanstd(nanmean(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){nMod}),2));
        end

        for x=1:numel(modalitiesInd{nInd})

            h=bar(x,dataMeans(x)); hold all
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

        nAov=findcolExact(aov.Effect,condNames{indEffect(nInd)});
        pMAIN=aov{nAov,6};

        amp=[max(yl)-min(yl)];
        isSignificant=0;
        pValues=ones(1,size(postHoc.(condNames{indEffect(nInd)}),1));
        pSelected=0;
        if pMAIN<pcritical(1) &  pSelected==0
            isSignificant=1;
            pValues=postHoc.(condNames{indEffect(nInd)}){:,5};
            pSelected=1;
        end

        if pMAIN<pcritical(end) &  pSelected==0
            pValues=postHoc.(condNames{indEffect(nInd)}){:,5};
            pSelected=1;
        end

        nColSignificant=ones(1,numel(dataMeans));
        nSignificant=1;
        if pSelected==1
            for i=1:numel(pValues)
                pV=pValues(i);
                if pV<=pcritical(end)
                    if yl(2)>0
                        if abs(dataMeans(order4ES.(condNames{indEffect(nInd)})(i,1)))>abs(dataMeans(order4ES.(condNames{indEffect(nInd)})(i,2)))
                            addPvalue(order4ES.(condNames{indEffect(nInd)})(i,1), dataMeans(order4ES.(condNames{indEffect(nInd)})(i,1))+0.05*nColSignificant(order4ES.(condNames{indEffect(nInd)})(i,1))*amp, pV, pcritical, colors{indEffect(nInd)}(order4ES.(condNames{indEffect(nInd)})(i,2),:))
                            nColSignificant(order4ES.(condNames{indEffect(nInd)})(i,1))=nColSignificant(order4ES.(condNames{indEffect(nInd)})(i,1))+1;
                        else
                            addPvalue(order4ES.(condNames{indEffect(nInd)})(i,2), dataMeans(order4ES.(condNames{indEffect(nInd)})(i,2))+0.05*nColSignificant(order4ES.(condNames{indEffect(nInd)})(i,2))*amp, pV, pcritical, colors{indEffect(nInd)}(order4ES.(condNames{indEffect(nInd)})(i,1),:))
                            nColSignificant(order4ES.(condNames{indEffect(nInd)})(i,2))=nColSignificant(order4ES.(condNames{indEffect(nInd)})(i,2))+1;

                        end
                    else
                        if abs(dataMeans(order4ES.(condNames{indEffect(nInd)})(i,1)))>abs(dataMeans(order4ES.(condNames{indEffect(nInd)})(i,2)))
                            addPvalue(order4ES.(condNames{indEffect(nInd)})(i,1), dataMeans(order4ES.(condNames{indEffect(nInd)})(i,1))-0.05*nColSignificant(order4ES.(condNames{indEffect(nInd)})(i,1))*amp, pV, pcritical, colors{indEffect(nInd)}(order4ES.(condNames{indEffect(nInd)})(i,2),:))
                            nColSignificant(order4ES.(condNames{indEffect(nInd)})(i,1))=nColSignificant(order4ES.(condNames{indEffect(nInd)})(i,1))+1;
                        else
                            addPvalue(order4ES.(condNames{indEffect(nInd)})(i,2), dataMeans(order4ES.(condNames{indEffect(nInd)})(i,2))-0.05*nColSignificant(order4ES.(condNames{indEffect(nInd)})(i,2))*amp, pV, pcritical, colors{indEffect(nInd)}(order4ES.(condNames{indEffect(nInd)})(i,1),:))
                            nColSignificant(order4ES.(condNames{indEffect(nInd)})(i,2))=nColSignificant(order4ES.(condNames{indEffect(nInd)})(i,2))+1;
                        end
                    end
                    nSignificant=nSignificant+1;
                end
            end
        end

        yl=yl*1.05.^max(nColSignificant);
        set(ax,'ylim',[min(yl) max(yl)])

        print('-dtiff',['-r' num2str(imageResolution)], fullfile(saveDir, 'text', condNames{indEffect(nInd)}, 'All RM'))
        close
        clear ax yl dataMeans dataSD

    end
end

% 1 IND 1RM
if numel(indEffect)>0 & numel(rmEffect)>0
    % Lines
    if plotLines
        for nRm=1:numel(effectRM)
            for nInd=1:numel(indEffect)

                f=figure('units','centimeters','position',[0 0 10+4*numel(modalitiesInd{nInd}) 5+9/16*4*numel(modalitiesInd{nInd})*numel(modalitiesRM{nRm})],'visible','off');

                for nMod=1:numel(allMod.(condNames{rmEffect(nRm)}))

                    subplot(numel(allMod.(condNames{rmEffect(nRm)})),1,nMod)

                    for x=1:numel(modalitiesInd{nInd})
                        dataMeans(x)=nanmean(nanmean(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){x})(:,col4means{nRm}(nMod,:)),2));
                        dataSD(x)=nanstd(nanmean(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){x})(:,col4means{nRm}(nMod,:)),2));
                    end

                    for x=1:numel(modalitiesInd{nInd})

                        h=bar(x,dataMeans(x)); hold all
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
                    if nMod==numel(allMod.(condNames{rmEffect(nRm)}))
                        xlabel(condNames{indEffect(nInd)})
                    end
                    ylabel(units)
                    yl(:,nMod)=ylim;
                    box off
                    ax{nMod}=gca;
                    ax{nMod}.XGrid='off';
                    ax{nMod}.YGrid='on';
                    title(allMod.(condNames{rmEffect(nRm)}){nMod})
                end

                nAovInt=findcolExact(aov.Effect,[condNames{indEffect(nInd)} ':' condNames{rmEffect(nRm)}]);
                nAov=findcolExact(aov.Effect,condNames{indEffect(nInd)});
                pINT=aov{nAovInt,6};
                pMAIN=aov{nAov,6};

                amp=[max(max(yl))-min(min(yl))];

                for nMod=1:numel(allMod.(condNames{rmEffect(nRm)}))
                    isSignificant=0;
                    set(f,'CurrentAxes',ax{nMod});

                    for x=1:numel(modalitiesInd{nInd})
                        dataMeans(x)=nanmean(nanmean(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){x})(:,col4means{nRm}(nMod,:)),2));
                        dataSD(x)=nanstd(nanmean(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){x})(:,col4means{nRm}(nMod,:)),2));
                    end

                    pValues=ones(1,size(postHoc.(condNames{indEffect(nInd)}),1));
                    pSelected=0;
                    if pINT<pcritical(1)
                        isSignificant=1;
                        rows=findcolExact(postHoc.([condNames{indEffect(nInd)} 'By' condNames{rmEffect(nRm)}]){:,1},allMod.(condNames{rmEffect(nRm)}){nMod});
                        pValues=postHoc.([condNames{indEffect(nInd)} 'By' condNames{rmEffect(nRm)}]){rows,6};
                        pSelected=1;
                    end
                    %                     if pMAIN<pcritical(1) &  pSelected==0
                    %                         isSignificant=1;
                    %                         pValues=postHoc.(condNames{indEffect(nInd)}){:,5};
                    %                         pSelected=1;
                    %                     end
                    if pINT<pcritical(end) &  pSelected==0
                        rows=findcolExact(postHoc.([condNames{indEffect(nInd)} 'By' condNames{rmEffect(nRm)}]){:,1},allMod.(condNames{rmEffect(nRm)}){nMod});
                        pValues=postHoc.([condNames{indEffect(nInd)} 'By' condNames{rmEffect(nRm)}]){rows,6};
                        pSelected=1;
                    end
                    %                     if pMAIN<pcritical(end) &  pSelected==0
                    %                         pValues=postHoc.(condNames{indEffect(nInd)}){:,5};
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
                                if yl(2)>0
                                    if pV<pcritical(1) & isSignificant==1
                                        hline(yl(2)+0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNames{indEffect(nInd)})(i,:),'lineWidth',1.5);
                                    else
                                        hline(yl(2)+0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNames{indEffect(nInd)})(i,:),'lineWidth',1.5);
                                    end
                                else
                                    if pV<pcritical(1) & isSignificant==1
                                        hline(yl(1)-0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNames{indEffect(nInd)})(i,:),'lineWidth',1.5);
                                    else
                                        hline(yl(1)-0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNames{indEffect(nInd)})(i,:),'lineWidth',1.5);
                                    end
                                end
                                nSignificant=nSignificant+1;
                            end
                        end
                    end
                end

                yl=yl*1.05.^max(nSignificant);
                for nMod=1:numel(allMod.(condNames{rmEffect(nRm)}))
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

            f=figure('units','centimeters','position',[0 0 10+4*numel(modalitiesInd{nInd}) 5+9/16*4*numel(modalitiesInd{nInd})*numel(modalitiesRM{nRm})],'visible','off');

            for nMod=1:numel(allMod.(condNames{rmEffect(nRm)}))

                subplot(numel(allMod.(condNames{rmEffect(nRm)})),1,nMod)

                for x=1:numel(modalitiesInd{nInd})
                    dataMeans(x)=nanmean(nanmean(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){x})(:,col4means{nRm}(nMod,:)),2));
                    dataSD(x)=nanstd(nanmean(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){x})(:,col4means{nRm}(nMod,:)),2));
                end

                for x=1:numel(modalitiesInd{nInd})

                    h=bar(x,dataMeans(x)); hold all
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
                if nMod==numel(allMod.(condNames{rmEffect(nRm)}))
                    xlabel(condNames{indEffect(nInd)})
                end
                ylabel(units)
                yl(:,nMod)=ylim;
                box off
                ax{nMod}=gca;
                ax{nMod}.XGrid='off';
                ax{nMod}.YGrid='on';
                title(allMod.(condNames{rmEffect(nRm)}){nMod})
            end

            nAovInt=findcolExact(aov.Effect,[condNames{indEffect(nInd)} ':' condNames{rmEffect(nRm)}]);
            nAov=findcolExact(aov.Effect,condNames{indEffect(nInd)});
            pINT=aov{nAovInt,6};
            pMAIN=aov{nAov,6};

            amp=[max(max(yl))-min(min(yl))];

            for nMod=1:numel(allMod.(condNames{rmEffect(nRm)}))

                isSignificant=0;
                set(f,'CurrentAxes',ax{nMod});

                for x=1:numel(modalitiesInd{nInd})
                    dataMeans(x)=nanmean(nanmean(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){x})(:,col4means{nRm}(nMod,:)),2));
                    dataSD(x)=nanstd(nanmean(data4plot.(condNames{indEffect(nInd)}).(allMod.(condNames{indEffect(nInd)}){x})(:,col4means{nRm}(nMod,:)),2));
                end

                pValues=ones(1,size(postHoc.(condNames{indEffect(nInd)}),1));
                pSelected=0;
                if pINT<pcritical(1)
                    isSignificant=1;
                    rows=findcolExact(postHoc.([condNames{indEffect(nInd)} 'By' condNames{rmEffect(nRm)}]){:,1},allMod.(condNames{rmEffect(nRm)}){nMod});
                    pValues=postHoc.([condNames{indEffect(nInd)} 'By' condNames{rmEffect(nRm)}]){rows,6};
                    pSelected=1;
                end
                %                 if pMAIN<pcritical(1) &  pSelected==0
                %                     isSignificant=1;
                %                     pValues=postHoc.(condNames{indEffect(nInd)}){:,5};
                %                     pSelected=1;
                %                 end
                if pINT<pcritical(end) &  pSelected==0
                    rows=findcolExact(postHoc.([condNames{indEffect(nInd)} 'By' condNames{rmEffect(nRm)}]){:,1},allMod.(condNames{rmEffect(nRm)}){nMod});
                    pValues=postHoc.([condNames{indEffect(nInd)} 'By' condNames{rmEffect(nRm)}]){rows,6};
                    pSelected=1;
                end
                %                 if pMAIN<pcritical(end) &  pSelected==0
                %                     pValues=postHoc.(condNames{indEffect(nInd)}){:,5};
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
                            if yl(2)>0
                                if abs(dataMeans(order4ES.(condNames{indEffect(nInd)})(i,1)))>abs(dataMeans(order4ES.(condNames{indEffect(nInd)})(i,2)))
                                    addPvalue(order4ES.(condNames{indEffect(nInd)})(i,1), dataMeans(order4ES.(condNames{indEffect(nInd)})(i,1))+0.05*nColSignificant(order4ES.(condNames{indEffect(nInd)})(i,1))*amp, pV, pcritical, colors{indEffect(nInd)}(order4ES.(condNames{indEffect(nInd)})(i,2),:))
                                    nColSignificant(order4ES.(condNames{indEffect(nInd)})(i,1))=nColSignificant(order4ES.(condNames{indEffect(nInd)})(i,1))+1;
                                else
                                    addPvalue(order4ES.(condNames{indEffect(nInd)})(i,2), dataMeans(order4ES.(condNames{indEffect(nInd)})(i,2))+0.05*nColSignificant(order4ES.(condNames{indEffect(nInd)})(i,2))*amp, pV, pcritical, colors{indEffect(nInd)}(order4ES.(condNames{indEffect(nInd)})(i,1),:))
                                    nColSignificant(order4ES.(condNames{indEffect(nInd)})(i,2))=nColSignificant(order4ES.(condNames{indEffect(nInd)})(i,2))+1;

                                end
                            else
                                if abs(dataMeans(order4ES.(condNames{indEffect(nInd)})(i,1)))>abs(dataMeans(order4ES.(condNames{indEffect(nInd)})(i,2)))
                                    addPvalue(order4ES.(condNames{indEffect(nInd)})(i,1), dataMeans(order4ES.(condNames{indEffect(nInd)})(i,1))-0.05*nColSignificant(order4ES.(condNames{indEffect(nInd)})(i,1))*amp, pV, pcritical, colors{indEffect(nInd)}(order4ES.(condNames{indEffect(nInd)})(i,2),:))
                                    nColSignificant(order4ES.(condNames{indEffect(nInd)})(i,1))=nColSignificant(order4ES.(condNames{indEffect(nInd)})(i,1))+1;
                                else
                                    addPvalue(order4ES.(condNames{indEffect(nInd)})(i,2), dataMeans(order4ES.(condNames{indEffect(nInd)})(i,2))-0.05*nColSignificant(order4ES.(condNames{indEffect(nInd)})(i,2))*amp, pV, pcritical, colors{indEffect(nInd)}(order4ES.(condNames{indEffect(nInd)})(i,1),:))
                                    nColSignificant(order4ES.(condNames{indEffect(nInd)})(i,2))=nColSignificant(order4ES.(condNames{indEffect(nInd)})(i,2))+1;
                                end
                            end
                            nSignificant=nSignificant+1;
                        end
                    end
                end
            end

            yl=yl*1.05.^max(nColSignificant);
            for nMod=1:numel(allMod.(condNames{rmEffect(nRm)}))
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
    if plotLines
        for nInd1=1:numel(indEffect)
            for nInd2=1:numel(indEffect)
                if nInd1~=nInd2

                    f=figure('units','centimeters','position',[0 0 10+4*numel(modalitiesInd{nInd1}) 5+9/16*4*numel(modalitiesInd{nInd1})*numel(modalitiesInd{nInd2})],'visible','off');

                    for nMod2=1:numel(modalitiesInd{nInd2})

                        subplot(numel(modalitiesInd{nInd2}),1,nMod2)

                        for nMod1=1:numel(modalitiesInd{nInd1})
                            dataMeans(nMod1)=nanmean(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNames{indEffect(nInd1)}){nMod1}).(allMod.(condNames{indEffect(nInd2)}){nMod2}),2));
                            dataSD(nMod1)=nanstd(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNames{indEffect(nInd1)}){nMod1}).(allMod.(condNames{indEffect(nInd2)}){nMod2}),2));
                        end

                        for x=1:numel(modalitiesInd{nInd1})

                            h=bar(x,dataMeans(x)); hold all
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
                        title(allMod.(condNames{indEffect(nInd2)}){nMod2})

                    end

                    nAovInt=findcolExact(aov.Effect,[condNames{indEffect(nInd1)} ':' condNames{indEffect(nInd2)}]);
                    if isempty(nAovInt)
                        nAovInt=findcolExact(aov.Effect,[condNames{indEffect(nInd2)} ':' condNames{indEffect(nInd1)}]);
                    end
                    nAov=findcolExact(aov.Effect,condNames{indEffect(nInd1)});
                    pINT=aov{nAovInt,6};
                    pMAIN=aov{nAov,6};

                    amp=[max(max(yl))-min(min(yl))];

                    for nMod2=1:numel(modalitiesInd{nInd2})

                        for nMod1=1:numel(modalitiesInd{nInd1})
                            dataMeans(nMod1)=nanmean(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNames{indEffect(nInd1)}){nMod1}).(allMod.(condNames{indEffect(nInd2)}){nMod2}),2));
                            dataSD(nMod1)=nanstd(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNames{indEffect(nInd1)}){nMod1}).(allMod.(condNames{indEffect(nInd2)}){nMod2}),2));
                        end

                        isSignificant=0;
                        set(f,'CurrentAxes',ax{nMod2});
                        pValues=ones(1,size(postHoc.(condNames{indEffect(nInd1)}),1));
                        pSelected=0;
                        if pINT<pcritical(1)
                            isSignificant=1;
                            phCut=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,1}, allMod.(condNames{indEffect(nInd2)}){nMod2});
                            pValues=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){phCut,6};
                            pSelected=1;
                        end
                        %                         if pMAIN<pcritical(1) &  pSelected==0
                        %                             isSignificant=1;
                        %                             pValues=postHoc.(condNames{indEffect(nInd1)}){:,5};
                        %                             pSelected=1;
                        %                         end
                        if pINT<pcritical(end) &  pSelected==0
                            phCut=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,1}, allMod.(condNames{indEffect(nInd2)}){nMod2});
                            pValues=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){phCut,6};
                            pSelected=1;
                        end
                        %                         if pMAIN<pcritical(end) &  pSelected==0
                        %                             pValues=postHoc.(condNames{indEffect(nInd1)}){:,5};
                        %                             pSelected=1;
                        %                         end

                        if pSelected==1
                            nSignificant=1;
                            for i=1:numel(pValues)
                                pV=pValues(i);
                                if pV<=pcritical(end)
                                    if yl(2,nMod2)>0
                                        if pV<pcritical(1) & isSignificant==1
                                            hline(yl(2,nMod2)+0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNames{indEffect(nInd1)})(i,:),'lineWidth',1.5);
                                        else
                                            hline(yl(2,nMod2)+0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNames{indEffect(nInd1)})(i,:),'lineWidth',1.5);
                                        end
                                    else
                                        if pV<pcritical(1) & isSignificant==1
                                            hline(yl(1,nMod2)-0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNames{indEffect(nInd1)})(i,:),'lineWidth',1.5);
                                        else
                                            hline(yl(1,nMod2)-0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNames{indEffect(nInd1)})(i,:),'lineWidth',1.5);
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

                f=figure('units','centimeters','position',[0 0 10+4*numel(modalitiesInd{nInd1}) 5+9/16*4*numel(modalitiesInd{nInd1})*numel(modalitiesInd{nInd2})],'visible','off');

                for nMod2=1:numel(modalitiesInd{nInd2})

                    subplot(numel(modalitiesInd{nInd2}),1,nMod2)

                    for nMod1=1:numel(modalitiesInd{nInd1})
                        dataMeans(nMod1)=nanmean(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNames{indEffect(nInd1)}){nMod1}).(allMod.(condNames{indEffect(nInd2)}){nMod2}),2));
                        dataSD(nMod1)=nanstd(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNames{indEffect(nInd1)}){nMod1}).(allMod.(condNames{indEffect(nInd2)}){nMod2}),2));
                    end

                    for x=1:numel(modalitiesInd{nInd1})

                        h=bar(x,dataMeans(x)); hold all
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
                    title(allMod.(condNames{indEffect(nInd2)}){nMod2})

                end

                nAovInt=findcolExact(aov.Effect,[condNames{indEffect(nInd1)} ':' condNames{indEffect(nInd2)}]);
                if isempty(nAovInt)
                    nAovInt=findcolExact(aov.Effect,[condNames{indEffect(nInd2)} ':' condNames{indEffect(nInd1)}]);
                end
                nAov=findcolExact(aov.Effect,condNames{indEffect(nInd1)});
                pINT=aov{nAovInt,6};
                pMAIN=aov{nAov,6};

                amp=[max(max(yl))-min(min(yl))];

                for nMod2=1:numel(modalitiesInd{nInd2})

                    for nMod1=1:numel(modalitiesInd{nInd1})
                        dataMeans(nMod1)=nanmean(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNames{indEffect(nInd1)}){nMod1}).(allMod.(condNames{indEffect(nInd2)}){nMod2}),2));
                        dataSD(nMod1)=nanstd(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNames{indEffect(nInd1)}){nMod1}).(allMod.(condNames{indEffect(nInd2)}){nMod2}),2));
                    end

                    if plotSD==1
                        dataMeans=dataMeans+sign(dataMeans).*dataSD;
                    end

                    isSignificant=0;
                    set(f,'CurrentAxes',ax{nMod2});
                    pValues=ones(1,size(postHoc.(condNames{indEffect(nInd1)}),1));
                    pSelected=0;
                    if pINT<pcritical(1)
                        isSignificant=1;
                        phCut=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,1}, allMod.(condNames{indEffect(nInd2)}){nMod2});
                        pValues=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){phCut,6};
                        pSelected=1;
                    end
                    %                     if pMAIN<pcritical(1) &  pSelected==0
                    %                         isSignificant=1;
                    %                         pValues=postHoc.(condNames{indEffect(nInd1)}){:,5};
                    %                         pSelected=1;
                    %                     end
                    if pINT<pcritical(end) &  pSelected==0
                        phCut=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,1}, allMod.(condNames{indEffect(nInd2)}){nMod2});
                        pValues=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){phCut,6};
                        pSelected=1;
                    end
                    %                     if pMAIN<pcritical(end) &  pSelected==0
                    %                         pValues=postHoc.(condNames{indEffect(nInd1)}){:,5};
                    %                         pSelected=1;
                    %                     end

                    nColSignificant=ones(1,numel(dataMeans));
                    nSignificant=1;
                    if pSelected==1
                        for i=1:numel(pValues)
                            pV=pValues(i);
                            if pV<=pcritical(end)
                                if yl(2)>0
                                    if abs(dataMeans(order4ES.(condNames{indEffect(nInd1)})(i,1)))>abs(dataMeans(order4ES.(condNames{indEffect(nInd1)})(i,2)))
                                        addPvalue(order4ES.(condNames{indEffect(nInd1)})(i,1), dataMeans(order4ES.(condNames{indEffect(nInd1)})(i,1))+0.05*nColSignificant(order4ES.(condNames{indEffect(nInd1)})(i,1))*amp, pV, pcritical, colors{indEffect(nInd1)}(order4ES.(condNames{indEffect(nInd1)})(i,2),:))
                                        nColSignificant(order4ES.(condNames{indEffect(nInd1)})(i,1))=nColSignificant(order4ES.(condNames{indEffect(nInd1)})(i,1))+1;
                                    else
                                        addPvalue(order4ES.(condNames{indEffect(nInd1)})(i,2), dataMeans(order4ES.(condNames{indEffect(nInd1)})(i,2))+0.05*nColSignificant(order4ES.(condNames{indEffect(nInd1)})(i,2))*amp, pV, pcritical, colors{indEffect(nInd1)}(order4ES.(condNames{indEffect(nInd1)})(i,1),:))
                                        nColSignificant(order4ES.(condNames{indEffect(nInd1)})(i,2))=nColSignificant(order4ES.(condNames{indEffect(nInd1)})(i,2))+1;
                                    end
                                else
                                    if abs(dataMeans(order4ES.(condNames{indEffect(nInd1)})(i,1)))>abs(dataMeans(order4ES.(condNames{indEffect(nInd1)})(i,2)))
                                        addPvalue(order4ES.(condNames{indEffect(nInd1)})(i,1), dataMeans(order4ES.(condNames{indEffect(nInd1)})(i,1))-0.05*nColSignificant(order4ES.(condNames{indEffect(nInd1)})(i,1))*amp, pV, pcritical, colors{indEffect(nInd1)}(order4ES.(condNames{indEffect(nInd1)})(i,2),:))
                                        nColSignificant(order4ES.(condNames{indEffect(nInd1)})(i,1))=nColSignificant(order4ES.(condNames{indEffect(nInd1)})(i,1))+1;
                                    else
                                        addPvalue(order4ES.(condNames{indEffect(nInd1)})(i,2), dataMeans(order4ES.(condNames{indEffect(nInd1)})(i,2))-0.05*nColSignificant(order4ES.(condNames{indEffect(nInd1)})(i,2))*amp, pV, pcritical, colors{indEffect(nInd1)}(order4ES.(condNames{indEffect(nInd1)})(i,1),:))
                                        nColSignificant(order4ES.(condNames{indEffect(nInd1)})(i,2))=nColSignificant(order4ES.(condNames{indEffect(nInd1)})(i,2))+1;
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
    if plotLines
        for nRm=1:numel(rmEffect)
            for nModRm=1:numel(modalitiesRM{nRm})
                for nInd1=1:numel(indEffect)
                    for nInd2=1:numel(indEffect)
                        if nInd1~=nInd2

                            f=figure('units','centimeters','position',[0 0 10+4*numel(modalitiesInd{nInd1}) 5+9/16*4*numel(modalitiesInd{nInd1})*numel(modalitiesInd{nInd2})],'visible','off');

                            for nMod2=1:numel(modalitiesInd{nInd2})

                                subplot(numel(modalitiesInd{nInd2}),1,nMod2)

                                for nMod1=1:numel(modalitiesInd{nInd1})
                                    dataMeans(nMod1)=nanmean(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNames{indEffect(nInd1)}){nMod1}).(allMod.(condNames{indEffect(nInd2)}){nMod2})(:,col4means{nRm}(nModRm,:)),2));
                                    dataSD(nMod1)=nanstd(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNames{indEffect(nInd1)}){nMod1}).(allMod.(condNames{indEffect(nInd2)}){nMod2})(:,col4means{nRm}(nModRm,:)),2));
                                end

                                for x=1:numel(modalitiesInd{nInd1})

                                    h=bar(x,dataMeans(x)); hold all
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
                                title(allMod.(condNames{indEffect(nInd2)}){nMod2})

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

                            nAov=findcolExact(aov.Effect,condNames{indEffect(nInd1)});
                            pINT3=aov{nAovInt3,6};
                            pINT2=aov{nAovInt2,6};
                            pINT=aov{nAovInt,6};
                            pMAIN=aov{nAov,6};

                            amp=[max(max(yl))-min(min(yl))];

                            for nMod2=1:numel(modalitiesInd{nInd2})

                                for nMod1=1:numel(modalitiesInd{nInd1})
                                    dataMeans(nMod1)=nanmean(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNames{indEffect(nInd1)}){nMod1}).(allMod.(condNames{indEffect(nInd2)}){nMod2})(:,col4means{nRm}(nModRm,:)),2));
                                    dataSD(nMod1)=nanstd(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNames{indEffect(nInd1)}){nMod1}).(allMod.(condNames{indEffect(nInd2)}){nMod2})(:,col4means{nRm}(nModRm,:)),2));
                                end

                                if plotSD==1
                                    dataMeans=dataMeans+sign(dataMeans).*dataSD;
                                end

                                isSignificant=0;
                                set(f,'CurrentAxes',ax{nMod2});
                                pValues=ones(1,size(postHoc.(condNames{indEffect(nInd1)}),1));
                                pSelected=0;
                                if pINT3<pcritical(1) & pSelected==0
                                    isSignificant=1;
                                    phCut1=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)} 'By' condNames{rmEffect(nRm)}]){:,1}, allMod.(condNames{rmEffect(nRm)}){nModRm});
                                    phCut2=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)} 'By' condNames{rmEffect(nRm)}]){:,2}, allMod.(condNames{indEffect(nInd2)}){nMod2});
                                    phCut=intersect(phCut1, phCut2);
                                    phCut=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)} 'By' condNames{rmEffect(nRm)}])(phCut,1:end);
                                    pValues=phCut{:,7};
                                    pSelected=1;
                                end
                                %                                 if any([pINT<pcritical(1) pINT2<pcritical(1)]) &  pSelected==0
                                %                                     isSignificant=1;
                                %                                     pSelected=1;
                                %                                     if pINT<pINT2
                                %                                         phCut=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{rmEffect(nRm)}]){:,1}, allMod.(condNames{rmEffect(nRm)}){nModRm});
                                %                                         pValues=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{rmEffect(nRm)}]){phCut,6};
                                %                                     else
                                %                                         phCut=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,1}, allMod.(condNames{indEffect(nInd2)}){nMod2});
                                %                                         pValues=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){phCut,6};
                                %                                     end
                                %                                 end
                                %                                 if pMAIN<pcritical(1) &  pSelected==0
                                %                                     isSignificant=1;
                                %                                     pValues=postHoc.(condNames{indEffect(nInd1)}){:,5};
                                %                                     pSelected=1;
                                %                                 end
                                if pINT3<pcritical(end) & pSelected==0
                                    isSignificant=1;
                                    phCut1=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)} 'By' condNames{rmEffect(nRm)}]){:,1}, allMod.(condNames{rmEffect(nRm)}){nModRm});
                                    phCut2=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)} 'By' condNames{rmEffect(nRm)}]){:,2}, allMod.(condNames{indEffect(nInd2)}){nMod2});
                                    phCut=intersect(phCut1, phCut2);
                                    phCut=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)} 'By' condNames{rmEffect(nRm)}])(phCut,1:end);
                                    pValues=phCut{:,7};
                                    pSelected=1;
                                end
                                %                                 if any([pINT<pcritical(end) pINT2<pcritical(end)]) &  pSelected==0
                                %                                     isSignificant=1;
                                %                                     pSelected=1;
                                %                                     if pINT<pINT2
                                %                                         phCut=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{rmEffect(nRm)}]){:,1}, allMod.(condNames{rmEffect(nRm)}){nModRm});
                                %                                         pValues=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{rmEffect(nRm)}]){phCut,6};
                                %                                     else
                                %                                         phCut=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,1}, allMod.(condNames{indEffect(nInd2)}){nMod2});
                                %                                         pValues=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){phCut,6};
                                %                                     end
                                %                                 end
                                %                                 if pMAIN<pcritical(end) &  pSelected==0
                                %                                     isSignificant=1;
                                %                                     pValues=postHoc.(condNames{indEffect(nInd1)}){:,5};
                                %                                     pSelected=1;
                                %                                 end

                                if pSelected==1
                                    nSignificant=1;
                                    for i=1:numel(pValues)
                                        pV=pValues(i);
                                        if pV<=pcritical(end)
                                            if yl(2,nMod2)>0
                                                if pV<pcritical(1) & isSignificant==1
                                                    hline(yl(2,nMod2)+0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNames{indEffect(nInd1)})(i,:),'lineWidth',1.5);
                                                else
                                                    hline(yl(2,nMod2)+0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNames{indEffect(nInd1)})(i,:),'lineWidth',1.5);
                                                end
                                            else
                                                if pV<pcritical(1) & isSignificant==1
                                                    hline(yl(1,nMod2)-0.05*nSignificant*amp,'linetype','-k','xLimits',order4ES.(condNames{indEffect(nInd1)})(i,:),'lineWidth',1.5);
                                                else
                                                    hline(yl(1,nMod2)-0.05*nSignificant*amp,'linetype','--k','xLimits',order4ES.(condNames{indEffect(nInd1)})(i,:),'lineWidth',1.5);
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

                        f=figure('units','centimeters','position',[0 0 10+4*numel(modalitiesInd{nInd1}) 5+9/16*4*numel(modalitiesInd{nInd1})*numel(modalitiesInd{nInd2})],'visible','off');

                        for nMod2=1:numel(modalitiesInd{nInd2})

                            subplot(numel(modalitiesInd{nInd2}),1,nMod2)

                            for nMod1=1:numel(modalitiesInd{nInd1})
                                dataMeans(nMod1)=nanmean(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNames{indEffect(nInd1)}){nMod1}).(allMod.(condNames{indEffect(nInd2)}){nMod2})(:,col4means{nRm}(nModRm,:)),2));
                                dataSD(nMod1)=nanstd(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNames{indEffect(nInd1)}){nMod1}).(allMod.(condNames{indEffect(nInd2)}){nMod2})(:,col4means{nRm}(nModRm,:)),2));
                            end

                            for x=1:numel(modalitiesInd{nInd1})

                                h=bar(x,dataMeans(x)); hold all
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
                            title(allMod.(condNames{indEffect(nInd2)}){nMod2})

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

                        nAov=findcolExact(aov.Effect,condNames{indEffect(nInd1)});
                        pINT3=aov{nAovInt3,6};
                        pINT2=aov{nAovInt2,6};
                        pINT=aov{nAovInt,6};
                        pMAIN=aov{nAov,6};

                        amp=[max(max(yl))-min(min(yl))];

                        for nMod2=1:numel(modalitiesInd{nInd2})

                            for nMod1=1:numel(modalitiesInd{nInd1})
                                dataMeans(nMod1)=nanmean(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNames{indEffect(nInd1)}){nMod1}).(allMod.(condNames{indEffect(nInd2)}){nMod2})(:,col4means{nRm}(nModRm,:)),2));
                                dataSD(nMod1)=nanstd(nanmean(data4plot.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]).(allMod.(condNames{indEffect(nInd1)}){nMod1}).(allMod.(condNames{indEffect(nInd2)}){nMod2})(:,col4means{nRm}(nModRm,:)),2));
                            end

                            if plotSD==1
                                dataMeans=dataMeans+sign(dataMeans).*dataSD;
                            end

                            isSignificant=0;
                            set(f,'CurrentAxes',ax{nMod2});
                            pValues=ones(1,size(postHoc.(condNames{indEffect(nInd1)}),1));
                            pSelected=0;
                            if pINT3<pcritical(1) & pSelected==0
                                isSignificant=1;
                                phCut1=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)} 'By' condNames{rmEffect(nRm)}]){:,1}, allMod.(condNames{rmEffect(nRm)}){nModRm});
                                phCut2=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)} 'By' condNames{rmEffect(nRm)}]){:,2}, allMod.(condNames{indEffect(nInd2)}){nMod2});
                                phCut=intersect(phCut1, phCut2);
                                phCut=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)} 'By' condNames{rmEffect(nRm)}])(phCut,1:end);
                                pValues=phCut{:,7};
                                pSelected=1;
                            end
                            %                             if any([pINT<pcritical(1) pINT2<pcritical(1)]) &  pSelected==0
                            %                                 isSignificant=1;
                            %                                 pSelected=1;
                            %                                 if pINT<pINT2
                            %                                     phCut=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{rmEffect(nRm)}]){:,1}, allMod.(condNames{rmEffect(nRm)}){nModRm});
                            %                                     pValues=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{rmEffect(nRm)}]){phCut,6};
                            %                                 else
                            %                                     phCut=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,1}, allMod.(condNames{indEffect(nInd2)}){nMod2});
                            %                                     pValues=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){phCut,6};
                            %                                 end
                            %                             end
                            %                             if pMAIN<pcritical(1) &  pSelected==0
                            %                                 isSignificant=1;
                            %                                 pValues=postHoc.(condNames{indEffect(nInd1)}){:,5};
                            %                                 pSelected=1;
                            %                             end
                            if pINT3<pcritical(end) & pSelected==0
                                isSignificant=1;
                                phCut1=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)} 'By' condNames{rmEffect(nRm)}]){:,1}, allMod.(condNames{rmEffect(nRm)}){nModRm});
                                phCut2=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)} 'By' condNames{rmEffect(nRm)}]){:,2}, allMod.(condNames{indEffect(nInd2)}){nMod2});
                                phCut=intersect(phCut1, phCut2);
                                phCut=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)} 'By' condNames{rmEffect(nRm)}])(phCut,1:end);
                                pValues=phCut{:,7};
                                pSelected=1;
                            end
                            %                             if any([pINT<pcritical(end) pINT2<pcritical(end)]) &  pSelected==0
                            %                                 isSignificant=1;
                            %                                 pSelected=1;
                            %                                 if pINT<pINT2
                            %                                     phCut=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{rmEffect(nRm)}]){:,1}, allMod.(condNames{rmEffect(nRm)}){nModRm});
                            %                                     pValues=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{rmEffect(nRm)}]){phCut,6};
                            %                                 else
                            %                                     phCut=findcol(postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){:,1}, allMod.(condNames{indEffect(nInd2)}){nMod2});
                            %                                     pValues=postHoc.([condNames{indEffect(nInd1)} 'By' condNames{indEffect(nInd2)}]){phCut,6};
                            %                                 end
                            %                             end
                            %                             if pMAIN<pcritical(end) &  pSelected==0
                            %                                 isSignificant=1;
                            %                                 pValues=postHoc.(condNames{indEffect(nInd1)}){:,5};
                            %                                 pSelected=1;
                            %                             end

                            nColSignificant=ones(1,numel(dataMeans));
                            nSignificant=1;
                            if pSelected==1
                                for i=1:numel(pValues)
                                    pV=pValues(i);
                                    if pV<=pcritical(end)
                                        if yl(2)>0
                                            if abs(dataMeans(order4ES.(condNames{indEffect(nInd1)})(i,1)))>abs(dataMeans(order4ES.(condNames{indEffect(nInd1)})(i,2)))
                                                addPvalue(order4ES.(condNames{indEffect(nInd1)})(i,1), dataMeans(order4ES.(condNames{indEffect(nInd1)})(i,1))+0.05*nColSignificant(order4ES.(condNames{indEffect(nInd1)})(i,1))*amp, pV, pcritical, colors{indEffect(nInd1)}(order4ES.(condNames{indEffect(nInd1)})(i,2),:))
                                                nColSignificant(order4ES.(condNames{indEffect(nInd1)})(i,1))=nColSignificant(order4ES.(condNames{indEffect(nInd1)})(i,1))+1;
                                            else
                                                addPvalue(order4ES.(condNames{indEffect(nInd1)})(i,2), dataMeans(order4ES.(condNames{indEffect(nInd1)})(i,2))+0.05*nColSignificant(order4ES.(condNames{indEffect(nInd1)})(i,2))*amp, pV, pcritical, colors{indEffect(nInd1)}(order4ES.(condNames{indEffect(nInd1)})(i,1),:))
                                                nColSignificant(order4ES.(condNames{indEffect(nInd1)})(i,2))=nColSignificant(order4ES.(condNames{indEffect(nInd1)})(i,2))+1;

                                            end
                                        else
                                            if abs(dataMeans(order4ES.(condNames{indEffect(nInd1)})(i,1)))>abs(dataMeans(order4ES.(condNames{indEffect(nInd1)})(i,2)))
                                                addPvalue(order4ES.(condNames{indEffect(nInd1)})(i,1), dataMeans(order4ES.(condNames{indEffect(nInd1)})(i,1))-0.05*nColSignificant(order4ES.(condNames{indEffect(nInd1)})(i,1))*amp, pV, pcritical, colors{indEffect(nInd1)}(order4ES.(condNames{indEffect(nInd1)})(i,2),:))
                                                nColSignificant(order4ES.(condNames{indEffect(nInd1)})(i,1))=nColSignificant(order4ES.(condNames{indEffect(nInd1)})(i,1))+1;
                                            else
                                                addPvalue(order4ES.(condNames{indEffect(nInd1)})(i,2), dataMeans(order4ES.(condNames{indEffect(nInd1)})(i,2))-0.05*nColSignificant(order4ES.(condNames{indEffect(nInd1)})(i,2))*amp, pV, pcritical, colors{indEffect(nInd1)}(order4ES.(condNames{indEffect(nInd1)})(i,1),:))
                                                nColSignificant(order4ES.(condNames{indEffect(nInd1)})(i,2))=nColSignificant(order4ES.(condNames{indEffect(nInd1)})(i,2))+1;
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