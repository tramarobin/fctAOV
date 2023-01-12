function unpackStructBase(structure)

fn = fieldnames(structure);

for i = 1:numel(fn)
    fni = string(fn(i));
    field = structure.(fni);
    if (isstruct(field))
        unpackStructBase(field);
        continue;
    end

        assignin('base', fni, field);
end



end