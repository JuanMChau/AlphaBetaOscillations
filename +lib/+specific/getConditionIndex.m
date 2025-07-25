function conditionIndex = getConditionIndex(condition,trialData)
%GETCONDITIONINDEX Summary of this function goes here
%   Detailed explanation goes here

conditionIndex = -1 * ones(size(condition,1),1);

for j = 1:size(condition,1)
    for i = 1:length(trialData)
        test = sum(ismember(trialData{i}{1},condition(j,1))) & ...
        sum(ismember(trialData{i}{2},condition(j,2))) & ...
        sum(ismember(trialData{i}{3},condition(j,3)));
        if (test)
            conditionIndex(j) = i;
            break;
        end
    end
end

end

