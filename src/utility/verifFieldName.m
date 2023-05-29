function fieldNameVerifed=verifFieldName(fieldName)

nTC={'/' '\' ':' '.' '<' '>' '?' '*' '|' '"' ' '}; % non tolerated characters

for i=1:numel(nTC)
    fieldName=strrep(fieldName,nTC{i},'');
end

fieldNameVerifed=fieldName;


end