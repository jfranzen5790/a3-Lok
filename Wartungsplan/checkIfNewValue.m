function newValue = checkIfNewValue(newValue)

    global allInputValues multipleInputValuesCounter;
    %for 
    currentValue = newValue;

    if   ~isempty(allInputValues)
        if  any(cellfun(@(x) isequal(x, currentValue), allInputValues))
            multipleInputValuesCounter = multipleInputValuesCounter + 1;
        else
            allInputValues{end + 1} = currentValue; 
        end
    else
        allInputValues{end + 1} = currentValue; 
    end
    
end

