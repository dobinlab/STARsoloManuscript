function outT = funTexBoldFirstWord(inT)

for ii=1:length(inT)
    inSpace = strfind(inT{ii},' ');
    outT{ii} = ['{\bf ' inT{ii}(1:inSpace-1) '}' inT{ii}(inSpace:end)];
end