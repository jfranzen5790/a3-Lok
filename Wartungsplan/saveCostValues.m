function null = saveCostValues(cellIndexCost, gueltig)
% Fkt. zum Abspeichern aller bisher berechneten Kosten
    %disp(size(allCosts(1) + 1));
    global allCosts costList gueltigCounter ungueltigCounter iterationsTillNewValidValue iterationsTillNewValidValue_ListOfallValues longestRuntime;
    allCosts{end+1} = cellIndexCost;
    costList(end+1) = cellIndexCost{2};
    if gueltig == 1 % ungültig!
        ungueltigCounter = ungueltigCounter + 1;
        
        iterationsTillNewValidValue = iterationsTillNewValidValue + 1;
        iterationsTillNewValidValue_ListOfallValues(end) = iterationsTillNewValidValue;
        
    else % gültig!
        gueltigCounter = gueltigCounter + 1;
        iterationsTillNewValidValue = 0;
        iterationsTillNewValidValue_ListOfallValues(end+1) = iterationsTillNewValidValue;
    end
    
    if iterationsTillNewValidValue > longestRuntime
       longestRuntime = iterationsTillNewValidValue; 
    end
end