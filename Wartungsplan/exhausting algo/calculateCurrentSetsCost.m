function [bestCost] = calculateCurrentSetsCost(currentInput, bestCost)
%TEST Summary of this function goes here
%   Detailed explanation goes here
    %currentInput = randomMaintenanceCombinations(value, :);
    
    
    
%     currentInputIndexVersion;
%     disp(currentInput);
     currentCost = (cost(currentInput));
     if currentCost < bestCost
         bestCost = currentCost;
         %bestComb = currentInput;
         %bestInput = currentInput;
     end
     %disp(bestCost);
% 
 %bestCost = bestCost + currentInput;
end