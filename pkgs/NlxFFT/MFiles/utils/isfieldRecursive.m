function isFieldResult = isfieldRecursive(inStruct, fieldName)
% isFieldResult = isfieldRecursive(inStruct, fieldName)
% inStruct is the name of the structure or an array of structures to search
% fieldName is the name of the field for which the function searches
isFieldResult = false;
if isempty(inStruct)
    return
end
f = fieldnames(inStruct(1));
for i=1:length(f)
    if(strcmp(f{i},strtrim(fieldName)))
        isFieldResult = true;
        return;
    elseif isstruct(inStruct(1).(f{i}))
            isFieldResult = isfieldRecursive(inStruct(1).(f{i}), fieldName);
            if isFieldResult
                return;
            end
    end
end