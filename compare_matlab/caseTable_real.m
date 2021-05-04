if ~exist('genome_kb','var')
    genome_kb='standard_1';
end

prefix1 = 'count/';
if ~exist('opt_alevinFry','var')
    opt_alevinFry='';
end

alevinSubdir='/gpl_knee/';

yesMult='mult:Yes';
noMult='mult:No';

caseTable = table;

caseTable = [caseTable;  {'CellRanger 5.0.1'                            'CellRanger'                                    ''                      'matrix.mtx'        'barcodes.tsv'         0                                       [prefix1 'CellRanger_5.0.1/'  genome1 '/standard/default/' sampleType      sample1 run1 'Run1/outs/filtered_feature_bc_matrix/']}];

caseTable = [caseTable;  {'STAR fullSA mult:No CR4'                     '{\bf STAR} full-SA'                            ''                      'matrix.mtx'        ''                     0                           [prefix1 'STAR_2.7.9x/'       genome1 '/fullSA/'       opt_STAR  '/' sampleType   '/'   sample1 '/' run1 'Solo.out/Gene/raw/']}];

caseTable = [caseTable;  {'STAR sparseSA mult:No CR4'                   '{\bf STAR} sparse-SA-3'                        ''                      'matrix.mtx'        ''                     0                           [prefix1 'STAR_2.7.9x/'       genome1 '/sparseSA3/'    opt_STAR sampleType      sample1 run1 'Solo.out/Gene/raw/']}];

genome_kb='standard_1';
caseTable = [caseTable;  {'kallisto-bustools mult:No'               ['{\bf kallisto|bustools} ' noMult]             ''                      'cells_x_genes.mtx' 'cells_x_genes.barcodes.txt'    1        [prefix1 'kbpy_0.25.0/'                      genome1 '/' genome_kb '/default' sampleType            sample1 run1 'counts_unfiltered/']}];
caseTable = [caseTable;  {'kallisto-bustools mult:Yes'              ['{\bf kallisto|bustools} ' yesMult]            ''                      'cells_x_genes.mtx' 'cells_x_genes.barcodes.txt'    1        [prefix1 'kbpy_0.25.0/'                      genome1 '/' genome_kb '/mult' sampleType               sample1 run1 'counts_unfiltered/']}];

%caseTable = [caseTable;  {'alevin-fry sketch trivial'                   '{\bf alevin-fry} sketch {\it trivial}'         ''                      'quants_mat.mtx'    'quants_mat_rows.txt'              1       [prefix1 'salmon-alevin-fry_1.4.0_0.1.0/'    genome1 '/standard/sketch_rad' opt_alevinFry sampleType   sample1 run1 alevinSubdir '/quant_trivial/alevin/']}];
caseTable = [caseTable;  {'alevin-fry sketch cr-like'                   ['{\bf alevin-fry} sketch ' noMult]             ''                      'quants_mat.mtx'    'quants_mat_rows.txt'              1       [prefix1 'salmon-alevin-fry_1.4.0_0.1.0/'    genome1 '/standard/sketch_rad' opt_alevinFry sampleType   sample1 run1 alevinSubdir '/quant_cr-like/alevin/']}];
caseTable = [caseTable;  {'alevin-fry sketch cr-like-em'                ['{\bf alevin-fry} sketch ' yesMult]            ''                      'quants_mat.mtx'    'quants_mat_rows.txt'              1       [prefix1 'salmon-alevin-fry_1.4.0_0.1.0/'    genome1 '/standard/sketch_rad' opt_alevinFry sampleType   sample1 run1 alevinSubdir '/quant_cr-like-em/alevin/']}];

%caseTable = [caseTable;  {'alevin-fry sel-align trivial'                '{\bf alevin-fry} sel-align {\it trivial}'      ''                      'quants_mat.mtx'    'quants_mat_rows.txt'          1       [prefix1 'salmon-alevin-fry_1.4.0_0.1.0/'    genome1 '/standard/rad' opt_alevinFry sampleType          sample1 run1 alevinSubdir '/quant_trivial/alevin/']}];
caseTable = [caseTable;  {'alevin-fry sel-align cr-like'                ['{\bf alevin-fry} sel-align ' noMult]          ''                      'quants_mat.mtx'    'quants_mat_rows.txt'          1       [prefix1 'salmon-alevin-fry_1.4.0_0.1.0/'    genome1 '/standard/rad' opt_alevinFry sampleType          sample1 run1 alevinSubdir '/quant_cr-like/alevin/']}];
caseTable = [caseTable;  {'alevin-fry sel-align cr-like-em'             ['{\bf alevin-fry} sel-align ' yesMult]         ''                      'quants_mat.mtx'    'quants_mat_rows.txt'          1       [prefix1 'salmon-alevin-fry_1.4.0_0.1.0/'    genome1 '/standard/rad' opt_alevinFry sampleType          sample1 run1 alevinSubdir '/quant_cr-like-em/alevin/']}];

%caseTable = [caseTable;  {'alevin-fry partial-decoy trivial'            '{\bf alevin-fry} partial-decoy {\it trivial}'  ''                      'quants_mat.mtx'    'quants_mat_rows.txt'          1       [prefix1 'salmon-alevin-fry_1.4.0_0.1.0/'    genome1 '/decoyPartial/rad' opt_alevinFry sampleType     sample1 run1 alevinSubdir '/quant_trivial/alevin/']}];
caseTable = [caseTable;  {'alevin-fry partial-decoy cr-like'            ['{\bf alevin-fry} partial-decoy ' noMult]      ''                      'quants_mat.mtx'    'quants_mat_rows.txt'          1       [prefix1 'salmon-alevin-fry_1.4.0_0.1.0/'    genome1 '/decoyPartial/rad' opt_alevinFry sampleType     sample1 run1 alevinSubdir '/quant_cr-like/alevin/']}];
caseTable = [caseTable;  {'alevin-fry partial-decoy cr-like-em'         ['{\bf alevin-fry} partial-decoy ' yesMult]     ''                      'quants_mat.mtx'    'quants_mat_rows.txt'          1       [prefix1 'salmon-alevin-fry_1.4.0_0.1.0/'    genome1 '/decoyPartial/rad' opt_alevinFry sampleType     sample1 run1 alevinSubdir '/quant_cr-like-em/alevin/']}];

%caseTable = [caseTable;  {'alevin-fry full-decoy trivial'               '{\bf alevin-fry} full-decoy {\it trivial}'     ''                      'quants_mat.mtx'    'quants_mat_rows.txt'          1       [prefix1 'salmon-alevin-fry_1.4.0_0.1.0/'    genome1 '/decoyFull/rad' opt_alevinFry sampleType        sample1 run1 alevinSubdir '/quant_trivial/alevin/']}];
caseTable = [caseTable;  {'alevin-fry full-decoy cr-like'               ['{\bf alevin-fry} full-decoy ' noMult]         ''                      'quants_mat.mtx'    'quants_mat_rows.txt'          1       [prefix1 'salmon-alevin-fry_1.4.0_0.1.0/'    genome1 '/decoyFull/rad' opt_alevinFry sampleType        sample1 run1 alevinSubdir '/quant_cr-like/alevin/']}];
caseTable = [caseTable;  {'alevin-fry full-decoy cr-like-em'            ['{\bf alevin-fry} full-decoy ' yesMult]        ''                      'quants_mat.mtx'    'quants_mat_rows.txt'          1       [prefix1 'salmon-alevin-fry_1.4.0_0.1.0/'    genome1 '/decoyFull/rad' opt_alevinFry sampleType        sample1 run1 alevinSubdir '/quant_cr-like-em/alevin/']}];

caseTable.Properties.VariableNames = {'tool'                            'label'                                         'geneFile'              'matFile'           'cbFile'                       'rowcol'  'dir' };

for ii=1:size(caseTable,1)
    caseTable.mat{ii} = [caseTable.dir{ii}   caseTable.matFile{ii}];
    if ~isempty(caseTable.cbFile{ii})
        caseTable.cb{ii} =  [caseTable.dir{ii}   caseTable.cbFile{ii}];
    else
        caseTable.cb{ii}='';
    end
    
    if startsWith(caseTable.tool{ii}, 'STAR')
        caseTable.lineStyle{ii}='-';
    elseif startsWith(caseTable.tool{ii}, 'kb_python')
        caseTable.lineStyle{ii}='--';    
    elseif startsWith(caseTable.tool{ii}, 'alevin-fry')
        caseTable.lineStyle{ii}=':';
    end

end

%
caseTable.labelSimple=strrep(caseTable.label,       '{\bf ', '');
caseTable.labelSimple=strrep(caseTable.labelSimple, '{\it ', '');
caseTable.labelSimple=strrep(caseTable.labelSimple, '}', '');
caseTable.labelSimple=strrep(caseTable.labelSimple, '\', '');
%caseTable.labelSimple=strrep(caseTable.labelSimple, ':', '');