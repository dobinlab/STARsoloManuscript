function toolColors = funToolsColors_Main(toolNames)

nTools=7;
vColor1=linspecer(nTools, 'qualitative');

for ii=1:length(toolNames)
    
    if startsWith(toolNames{ii}, 'STAR') && contains(toolNames{ii}, 'parse')
        toolColors{ii}=vColor1(1,:);
    elseif startsWith(toolNames{ii}, 'STAR')
        toolColors{ii}=vColor1(2,:);

    elseif startsWith(toolNames{ii}, 'alevin-fry sketch') || startsWith(toolNames{ii}, 'alevin-fry~sketch')
        toolColors{ii}=vColor1(4,:);
    elseif startsWith(toolNames{ii}, 'alevin-fry sel-align') || startsWith(toolNames{ii}, 'alevin-fry~sel-align')
        toolColors{ii}=vColor1(5,:);
    elseif startsWith(toolNames{ii}, 'alevin-fry partial-decoy') || startsWith(toolNames{ii}, 'alevin-fry~partial-decoy')
        toolColors{ii}=vColor1(6,:);
    elseif startsWith(toolNames{ii}, 'alevin-fry full-decoy') || startsWith(toolNames{ii}, 'alevin-fry~full-decoy')
        toolColors{ii}=vColor1(7,:);       
    elseif startsWith(toolNames{ii}, 'kb') || startsWith(toolNames{ii}, 'kallisto')
        toolColors{ii}=vColor1(3,:);  
        
    elseif startsWith(toolNames{ii}, 'CellRanger')
        toolColors{ii}=[0 0 0];                
    end

end