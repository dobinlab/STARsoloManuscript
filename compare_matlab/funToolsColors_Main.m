function [toolColors, toolMarker] = funToolsColors_Main(toolNames)

nTools=7;
vColor1=linspecer(nTools, 'qualitative');

for ii=1:length(toolNames)
    
    if startsWith(toolNames{ii}, 'STAR') && ~contains(toolNames{ii}, 'parse')
        toolColors{ii}=[255,48,39]/255;         % red
        toolMarker{ii}='p';
    elseif startsWith(toolNames{ii}, 'STAR')
        toolColors{ii}=[252,141,89]/255;     % light red
        toolMarker{ii}='*';
    %elseif startsWith(toolNames{ii}, 'alevin-fry sketch') || startsWith(toolNames{ii}, 'alevin-fry~sketch')
    elseif contains(toolNames{ii},'sketch')
        toolColors{ii}=[90,180,172]/255;
        toolMarker{ii}='>';
    %elseif startsWith(toolNames{ii}, 'alevin-fry sel-align') || startsWith(toolNames{ii}, 'alevin-fry~sel-align')
    elseif contains(toolNames{ii},'sel-align')
        toolColors{ii}=[153,142,195]/255;
        toolMarker{ii}='<';        
    %elseif startsWith(toolNames{ii}, 'alevin-fry partial-decoy') || startsWith(toolNames{ii}, 'alevin-fry~partial-decoy')
    elseif contains(toolNames{ii},'partial-decoy')
        toolColors{ii}=[145,191,219]/255;  
        toolMarker{ii}='v';
    %elseif startsWith(toolNames{ii}, 'alevin-fry full-decoy') || startsWith(toolNames{ii}, 'alevin-fry~full-decoy')
    elseif contains(toolNames{ii},'full-decoy')
        toolColors{ii}=[69,117,180]/255;        %blue
        toolMarker{ii}='^';
    elseif startsWith(toolNames{ii}, 'kb') || contains(toolNames{ii}, 'allisto')
        toolColors{ii}=[27,120,55]/255;  % green
        toolMarker{ii}='o';
    elseif startsWith(toolNames{ii}, 'CellRanger')
        toolColors{ii}=[0 0 0];
        toolMarker{ii}='s';
    end

end