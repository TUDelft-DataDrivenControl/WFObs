function [ Sout ] = mergeStruct( S1,S2 )
    Sout = S1; % Copy values from S1
    fieldNames = fieldnames(S2);
    for j = 1:length(fieldNames)
        fN = (fieldNames{j});
        if strcmp(class(S2.(fN)), 'struct')
            Sout.(fN) = mergeStruct(Sout.(fN),S2.(fN)); % Go into substruct
        else
            Sout.(fN) = S2.(fN); % Overwrite value
        end
    end
end