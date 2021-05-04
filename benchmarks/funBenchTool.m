function t=funBenchTool(name1,genome1,opt1)

if startsWith(name1,'CellRanger')
    t='CellRanger';
    return;
elseif startsWith(name1,'kbpy')
    t='Kallisto';
    return;
elseif startsWith(name1,'STAR')    
    t='STARsolo';
    if strcmp(genome1,'fullSA')
        t=[t '_fullSA'];
    elseif strcmp(genome1,'sparseSA3')
        t=[t '_sparseSA'];
    end
        
elseif startsWith(name1,'salmon')
    t='Alevin';
    if strcmp(genome1,'decoyFull')
        t=[t '_full-decoy'];
    elseif strcmp(genome1,'decoyPartial')
        t=[t '_partial-decoy'];
    elseif strcmp(genome1,'standard') && strcmp(opt1,'rad')
        t=[t '_sel-align'];
    elseif strcmp(genome1,'standard') && strcmp(opt1,'sketch_rad')
        t=[t '_sketch'];
    end
    
else
    disp(name1);
end
