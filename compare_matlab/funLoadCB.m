function cb1=funLoadCB(fileIn)

fin1=fopen(fileIn);
cb1=textscan(fin1,'%s %*[^\n]');
cb1=vertcat(cb1{1}{:});
fclose(fin1);

if cb1(1,end-1)=='-'
    cb1=cb1(:,1:end-2);
end