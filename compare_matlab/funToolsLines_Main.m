function toolLines = funToolsLines_Main(toolNames)

for ii=1:length(toolNames)
    
    if contains(toolNames{ii}, 'STAR')
        toolLines{ii}='-';
    elseif startsWith(toolNames{ii}, 'kb') || contains(toolNames{ii}, 'allisto')
        toolLines{ii}='-.';    
    elseif contains(toolNames{ii}, 'levin')
        if contains(toolNames{ii}, 'decoy')
            toolLines{ii}='--';
        else
            toolLines{ii}=':';
        end
    elseif contains(toolNames{ii}, 'CellRanger')
        toolLines{ii}=':';
    end

end