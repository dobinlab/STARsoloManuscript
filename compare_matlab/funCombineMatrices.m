function [Ma,cb]=funCombineMatrices(M,Ma,cb)

%% find interesection of CBs
for ii=1:length(M)
    if ~isempty(M{ii})
        cb=cb & sum(M{ii})>0;
    end
end

disp(['Number of intersection cells:' num2str(nnz(cb))]);

for ii=1:length(M)
    fprintf(1, '%i ', ii);
    if ~isempty(M{ii}) && ( ii>length(Ma) || isempty(Ma{ii}) )
        Ma{ii}=M{ii}(:,cb);
    end
end
