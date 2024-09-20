function tableStats=createTableStats(tableStats, format, stats, tables, varName)

if numel(stats.condNames)==2

    e1=string(stats.condNames{1});
    e2=string(stats.condNames{2});
    c1=string(stats.cond4effect{1});
    c2=string(stats.cond4effect{2});
    nParticipants=numel(stats.ID);

    indEffect=any(find(~stats.isRM));
    if indEffect
        c1=unique(c1,'stable')';
    end


    if isempty(tableStats)
        tableStats=array2table({e2, e1});
        switch indEffect
            case 0
                tableStats=[tableStats; table([repmat(c2{1},numel(c1)+1,1); repmat(c2{2},numel(c1)+1,1); repmat("ALL",numel(c1),1)], [repmat([string(c1'); "ALL"],numel(c2),1); string(c1')])];
            case 1
                tableStats=[tableStats; table([repmat(c2{1},numel(c1)+1,1); repmat(c2{2},numel(c1)+1,1); repmat("ALL",numel(c1)+1,1)], [repmat([string(c1'); "ALL"],numel(c2)+1,1)])];
        end
        tableStats=[tableStats; table([e1; e2; e1 + " x " + e2], nan(3,1))];
    end

    loop=2;
    vt={varName};

    switch indEffect

        case 0

            for m2=1:numel(c2)+1
                for m1=1:numel(c1)+1

                    if m2<numel(c2)+1 & m1<numel(c1)+1

                        vt{loop,1}=sprintf([format ' \x00B1 ' format], tables.means.allData{nParticipants+2,[c1{m1} ' ' c2{m2}]}, tables.means.allData{nParticipants+3,[c1{m1} ' ' c2{m2}]});
                        loop=loop+1;

                    elseif m1>numel(c1) & m2<numel(c2)+1

                        vt{loop,1}=sprintf([format ' \x00B1 ' format], tables.means.(verifFieldName(e2)){nParticipants+2,[c2{m2}]}, tables.means.(verifFieldName(e2)){nParticipants+3,[c2{m2}]});
                        loop=loop+1;

                    elseif m2>numel(c2) & m1<numel(c1)+1

                        vt{loop,1}=sprintf([format ' \x00B1 ' format], tables.means.(verifFieldName(e1)){nParticipants+2,[c1{m1}]}, tables.means.(verifFieldName(e1)){nParticipants+3,[c1{m1}]});
                        loop=loop+1;

                    end

                end
            end

        case 1

            for m2=1:numel(stats.cond4effect{2})+1

                if m2<=numel(stats.cond4effect{2})

                    vt{loop,1}=sprintf([format ' \x00B1 ' format], tables.means.allData{nParticipants+5,c2{m2}}, tables.means.allData{nParticipants+6,c2{m2}});
                    vt{loop+1,1}=sprintf([format ' \x00B1 ' format], tables.means.allData{nParticipants+8,c2{m2}}, tables.means.allData{nParticipants+9,c2{m2}});
                    vt{loop+2,1}=sprintf([format ' \x00B1 ' format], tables.means.allData{nParticipants+2,c2{m2}}, tables.means.allData{nParticipants+3,c2{m2}});
                    loop=loop+3;

                else

                    vt{loop,1}=sprintf([format ' \x00B1 ' format], tables.means.allData{nParticipants+5,'ALL RM'}, tables.means.allData{nParticipants+6,'ALL RM'});
                    vt{loop+1,1}=sprintf([format ' \x00B1 ' format], tables.means.allData{nParticipants+8,'ALL RM'}, tables.means.allData{nParticipants+9,'ALL RM'});
                    vt{loop+2,1}=sprintf([format ' \x00B1 ' format], tables.means.allData{nParticipants+2,'ALL RM'}, tables.means.allData{nParticipants+3,'ALL RM'});
                    loop=loop+3;

                end
            end
    end

    for e=1:3
        if tables.aov{e,6}>0.001
            vt{loop,1}=sprintf('F=%.3g, p=%.3f', tables.aov{e,5}, tables.aov{e,6});
        else
            vt{loop,1}=sprintf('F=%.3g, p<%.3f', tables.aov{e,5}, 0.001);
        end
        vt{loop+1,1}=sprintf('%.4f', tables.aov{e,6});
        loop=loop+2;
    end

    tableStats(:,size(tableStats,2)+1)=vt;

elseif numel(stats.condNames)==1

    e1=string(stats.condNames{1});
    c1=string(stats.cond4effect{1});
    nParticipants=numel(stats.ID);

    indEffect=any(find(~stats.isRM));
    if indEffect
        c1=unique(c1,'stable')';
    end

    if isempty(tableStats)
        tableStats=array2table([e1; c1'; "ALL"; [e1 + " aov"]; [e1 + " p-value"]]);
    end

    loop=2;
    vt={varName};

    if indEffect

        for m1=1:numel(c1)+1

            if m1<numel(c1)+1

                vt{loop,1}=sprintf([format ' \x00B1 ' format], tables.means.allData{nParticipants+5+(m1-1)*3,3}, tables.means.allData{nParticipants+6+(m1-1)*3,3});
                loop=loop+1;

            else

                vt{loop,1}=sprintf([format ' \x00B1 ' format], tables.means.allData{nParticipants+2,3}, tables.means.allData{nParticipants+3,3});
                loop=loop+1;

            end

        end

    else

        for m1=1:numel(c1)+1

            vt{loop,1}=sprintf([format ' \x00B1 ' format], tables.means.allData{nParticipants+2,m1+1}, tables.means.allData{nParticipants+3,m1+1});
            loop=loop+1;

        end

    end

    if indEffect
        if tables.aov{1,6}{1}>0.001
            vt{loop,1}=sprintf('F=%.3g, p=%.3f', tables.aov{1,5}{1}, tables.aov{1,6}{1});
        else
            vt{loop,1}=sprintf('F=%.3g, p<%.3f', tables.aov{1,5}{1}, 0.001);
        end
        vt{loop+1,1}=sprintf('%.5f', tables.aov{1,6}{1});

    else

        if tables.aov{1,6}>0.001
            vt{loop,1}=sprintf('F=%.3g, p=%.3f', tables.aov{1,5}, tables.aov{1,6});
        else
            vt{loop,1}=sprintf('F=%.3g, p<%.3f', tables.aov{1,5}, 0.001);
        end
        vt{loop+1,1}=sprintf('%.5f', tables.aov{1,6});

    end

    tableStats(:,size(tableStats,2)+1)=vt;

end

end


