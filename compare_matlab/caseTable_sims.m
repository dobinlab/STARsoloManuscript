prefix1='count/';

yesMult='multYes';
noMult='multNo';

caseTable = table;
caseTable = [caseTable;  {'Simulated Truth'                 '{\bf Simulated Truth}'               'truth.mtx'          ''      0         ['samples' sampleType sample1 '/']}];

caseTable = [caseTable;  {'STARsolo fullSA mult:No ENCODE'          ['{\bf STARsolo} ' noMult ' ENCODE']                'matrix.mtx'                    ''      0         [prefix1 'STAR_2.7.8a/' genome1 '/fullSA/'       '10X_noSAM_sims_mult_ENCODE/' sampleType      sample1 run1 'Solo.out/Gene/raw/']}];
caseTable = [caseTable;  {'STARsolo fullSA mult:Uniform ENCODE'     ['{\bf STARsolo} ' yesMult ':Uniform' ' ENCODE']     'UniqueAndMult-Uniform.mtx'     ''      0         [prefix1 'STAR_2.7.8a/' genome1 '/fullSA/'       '10X_noSAM_sims_mult_ENCODE/' sampleType      sample1 run1 'Solo.out/Gene/raw/']}];
caseTable = [caseTable;  {'STARsolo fullSA mult:Rescue ENCODE'      ['{\bf STARsolo} ' yesMult ':Recsue' ' ENCODE']      'UniqueAndMult-Rescue.mtx'      ''      0         [prefix1 'STAR_2.7.8a/' genome1 '/fullSA/'       '10X_noSAM_sims_mult_ENCODE/' sampleType      sample1 run1 'Solo.out/Gene/raw/']}];
caseTable = [caseTable;  {'STARsolo fullSA mult:EM ENCODE'          ['{\bf STARsolo} ' yesMult ':EM' ' ENCODE']          'UniqueAndMult-EM.mtx'          ''      0         [prefix1 'STAR_2.7.8a/' genome1 '/fullSA/'       '10X_noSAM_sims_mult_ENCODE/' sampleType      sample1 run1 'Solo.out/Gene/raw/']}];


caseTable = [caseTable;  {'kallisto-bustools mult:No'                   ['{\bf kallisto|bustools} ' noMult]               'cells_x_genes.mtx'  'cells_x_genes.barcodes.txt'    1         [prefix1 'kbpy_0.25.0/'                      genome1 '/standard_1/default' sampleType            sample1 run1 'counts_unfiltered/']}];
caseTable = [caseTable;  {'kallisto-bustools mult:Yes'                  ['{\bf kallisto|bustools} ' yesMult]              'cells_x_genes.mtx'  'cells_x_genes.barcodes.txt'    1         [prefix1 'kbpy_0.25.0/'                      genome1 '/standard_1/mult' sampleType               sample1 run1 'counts_unfiltered/']}];

% caseTable = [caseTable;  {'alevin-fry sketch trivial'           '{\bf alevin-fry} sketch {\it trivial}'              'quants_mat.mtx' 'quants_mat_rows.txt' 1         [prefix1 'salmon-alevin-fry_1.4.0_0.1.0/'    genome1 '/standard/sketch_rad' opt_alevinFry sampleType   sample1 run1 'gpl_knee/quant_trivial/alevin/']}];
caseTable = [caseTable;  {'alevin-fry sketch cr-like'           ['{\bf alevin-fry} sketch ' noMult]              'quants_mat.mtx' 'quants_mat_rows.txt' 1         [prefix1 'salmon-alevin-fry_1.4.0_0.1.0/'    genome1 '/standard/sketch_rad' opt_alevinFry sampleType   sample1 run1 'gpl_knee/quant_cr-like/alevin/']}];
caseTable = [caseTable;  {'alevin-fry sketch cr-like-em'        ['{\bf alevin-fry} sketch ' yesMult]           'quants_mat.mtx' 'quants_mat_rows.txt' 1         [prefix1 'salmon-alevin-fry_1.4.0_0.1.0/'    genome1 '/standard/sketch_rad' opt_alevinFry sampleType   sample1 run1 'gpl_knee/quant_cr-like-em/alevin/']}];

% caseTable = [caseTable;  {'alevin-fry sel-align trivial'        '{\bf alevin-fry} sel-align {\it trivial}'           'quants_mat.mtx'     'quants_mat_rows.txt'           1         [prefix1 'salmon-alevin-fry_1.4.0_0.1.0/'    genome1 '/standard/rad' opt_alevinFry sampleType          sample1 run1 'gpl_knee/quant_trivial/alevin/']}];
caseTable = [caseTable;  {'alevin-fry sel-align cr-like'        ['{\bf alevin-fry} sel-align ' noMult]           'quants_mat.mtx'     'quants_mat_rows.txt'           1         [prefix1 'salmon-alevin-fry_1.4.0_0.1.0/'    genome1 '/standard/rad' opt_alevinFry sampleType          sample1 run1 'gpl_knee/quant_cr-like/alevin/']}];
caseTable = [caseTable;  {'alevin-fry sel-align cr-like-em'     ['{\bf alevin-fry} sel-align ' yesMult]        'quants_mat.mtx'     'quants_mat_rows.txt'           1         [prefix1 'salmon-alevin-fry_1.4.0_0.1.0/'    genome1 '/standard/rad' opt_alevinFry sampleType          sample1 run1 'gpl_knee/quant_cr-like-em/alevin/']}];

% caseTable = [caseTable;  {'alevin-fry partial-decoy trivial'    '{\bf alevin-fry} partial-decoy {\it trivial}'     'quants_mat.mtx'     'quants_mat_rows.txt'           1         [prefix1 'salmon-alevin-fry_1.4.0_0.1.0/'    genome1 '/decoyPartial/rad' opt_alevinFry sampleType     sample1 run1 'gpl_knee/quant_trivial/alevin/']}];
caseTable = [caseTable;  {'alevin-fry partial-decoy cr-like'    ['{\bf alevin-fry} partial-decoy ' noMult]     'quants_mat.mtx'     'quants_mat_rows.txt'           1         [prefix1 'salmon-alevin-fry_1.4.0_0.1.0/'    genome1 '/decoyPartial/rad' opt_alevinFry sampleType     sample1 run1 'gpl_knee/quant_cr-like/alevin/']}];
caseTable = [caseTable;  {'alevin-fry partial-decoy cr-like-em' ['{\bf alevin-fry} partial-decoy ' yesMult]  'quants_mat.mtx'     'quants_mat_rows.txt'           1         [prefix1 'salmon-alevin-fry_1.4.0_0.1.0/'    genome1 '/decoyPartial/rad' opt_alevinFry sampleType     sample1 run1 'gpl_knee/quant_cr-like-em/alevin/']}];

% caseTable = [caseTable;  {'alevin-fry full-decoy trivial'       '{\bf alevin-fry} full-decoy {\it trivial}'      'quants_mat.mtx'     'quants_mat_rows.txt'           1         [prefix1 'salmon-alevin-fry_1.4.0_0.1.0/'    genome1 '/decoyFull/rad' opt_alevinFry sampleType        sample1 run1 'gpl_knee/quant_trivial/alevin/']}];
caseTable = [caseTable;  {'alevin-fry full-decoy cr-like'       ['{\bf alevin-fry} full-decoy ' noMult]      'quants_mat.mtx'     'quants_mat_rows.txt'           1         [prefix1 'salmon-alevin-fry_1.4.0_0.1.0/'    genome1 '/decoyFull/rad' opt_alevinFry sampleType        sample1 run1 'gpl_knee/quant_cr-like/alevin/']}];
caseTable = [caseTable;  {'alevin-fry full-decoy cr-like-em'    ['{\bf alevin-fry} full-decoy ' yesMult]   'quants_mat.mtx'     'quants_mat_rows.txt'           1         [prefix1 'salmon-alevin-fry_1.4.0_0.1.0/'    genome1 '/decoyFull/rad' opt_alevinFry sampleType        sample1 run1 'gpl_knee/quant_cr-like-em/alevin/']}];

caseTable.Properties.VariableNames = {'tool'                     'label'                                'matFile'            'cbFile'                        'rowcol'      'dir'                                             };

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

caseTable.geneFile = repmat({''},size(caseTable,1),1);