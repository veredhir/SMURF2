function unpackStruct (filePath)
    if(isa(filePath,'char'))
        structure = load(filePath);
    else
        structure = filePath;
    end
    fn = fieldnames(structure);
    for i = 1:numel(fn)
        fni = string(fn(i));
        structure.(fni);
        field = structure.(fni);
        if (isstruct(field))
            unpackStruct(field);
            continue;
        end
        assignin('base', fni, field);
    end
end
