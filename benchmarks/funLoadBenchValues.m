function allValues=funLoadBenchValues(file1, toolNames, yesTimeConv)
dbstop if error;

wtT = readtable(file1);

%%
wtT.Properties.VariableNames={'name', 'genome', 'opt', 'threads', 'run', 'quant', 'values'};

%%
for ii=1:height(wtT)
    wtT.tool{ii} = funBenchTool(wtT.name{ii}, wtT.genome{ii}, wtT.opt{ii});
end

%%
if yesTimeConv
    for ii=1:height(wtT)
        vt = sscanf(wtT.values{ii},'%f:%f:%f');
        vt = [zeros(3-length(vt),1); vt];
        wtT.seconds(ii)=sum(vt.*[3600;60;1]);
    end
end

%%

runID={'b01', 'b02', 'b03', 'b04', 'b05' };

nThreads=4:4:20;
%
allValues = zeros(length(nThreads),length(toolNames),length(runID));

for itool=1:length(toolNames)
    fprintf(1, '%i ', itool);
    for ithr=1:length(nThreads)
        for ir=1:length(runID)
            vv=strcmp(toolNames{itool}, wtT.tool) & strcmp(runID{ir}, wtT.run) & nThreads(ithr)==wtT.threads;
            if yesTimeConv
                allValues(ithr, itool, ir)=sum(wtT.seconds(vv));
            elseif nnz(vv)>0
                allValues(ithr, itool, ir)=max(wtT.values(vv));
            end
        end
    end
end