function [RandomSelection] = RandomSelection(EdgeTabulation,RandomNode)
%Function to select a random node based on the tabulated probabilities

RandomSelection=0;
for j=1:numel(EdgeTabulation(:,2))
    
    if sum(EdgeTabulation(1:j,2))>=RandomNode
        RandomSelection=j;
    end
    
    %Break loop if random F has been found
    if RandomSelection>0
        break
    end
end

end