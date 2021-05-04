function funWriteMatrixTable(file1, sheet1, mat1, labelRow, labelCol)

C = cell(size(mat1)+1);
C(2:end,1) = labelRow;
C(1,2:end) = labelCol(1:size(mat1,2));

C(2:end,2:end) = arrayfun(@num2str, mat1, 'UniformOutput', 0);

writetable(cell2table(C), [file1 '.xlsx'], 'Sheet', sheet1, 'WriteVariableNames', false)