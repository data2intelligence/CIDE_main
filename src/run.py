#!/usr/bin/env python
import os, sys, pathlib, re, pandas, numpy, warnings
import statsmodels.api as sm
import hashlib
import CytoSig
import data_significance
import sklearn.metrics as metrics

from scipy import stats
from glob import glob
from statsmodels.stats.multitest import multipletests
from statsmodels.tools.sm_exceptions import MissingDataError

base_path = pathlib.Path(__file__).parent.absolute()
base_path = os.path.dirname(base_path)
data_path = os.path.join(base_path, 'data')
input_path = os.path.join(data_path, 'data_open')

lst_survival = []

def load_survival_list(data_type):
    # IO region
    IO_lst = [
        ['Pancreatic_Nivolumab_Padron2022', 'OS_Nivo+Chemo'],
        ['Pancreatic_Nivolumab_Padron2022', 'PFS_Nivo+Chemo'],
        
        ['Pancreatic_Nivolumab_Padron2022', 'OS_Nivo+Sotiga+Chemo'],
        ['Pancreatic_Nivolumab_Padron2022', 'PFS_Nivo+Sotiga+Chemo'],
        
        #['Pancreatic_Nivolumab_Padron2022', 'OS_Sotiga+Chemo'],
        #['Pancreatic_Nivolumab_Padron2022', 'PFS_Sotiga+Chemo'],
        
        ['Pancreatic_Nivolumab_Padron2022.Liver', 'OS_Nivo+Chemo'],
        ['Pancreatic_Nivolumab_Padron2022.Liver', 'PFS_Nivo+Chemo'],
        
        ['Pancreatic_Nivolumab_Padron2022.Liver', 'OS_Nivo+Sotiga+Chemo'],
        ['Pancreatic_Nivolumab_Padron2022.Liver', 'PFS_Nivo+Sotiga+Chemo'],
        
        ['Pancreatic_Nivolumab_Padron2022.Liver', 'OS_Sotiga+Chemo'],
        ['Pancreatic_Nivolumab_Padron2022.Liver', 'PFS_Sotiga+Chemo'],
        
        ['NSCLC_PD1orPDL1_Jung2019', 'PFS'],
        
        ['Melanoma_Ipilimumab_VanAllen2015', 'OS'],
        ['Melanoma_Ipilimumab_VanAllen2015', 'PFS'],
        
        ['Melanoma_PD1_Hugo2016', 'OS'],
        
        ['Melanoma_Nivolumab_Riaz2017', 'OS_Naive'],
        ['Melanoma_Nivolumab_Riaz2017', 'PFS_Naive'],
        
        ['Melanoma_Nivolumab_Riaz2017', 'OS_Prog'],
        ['Melanoma_Nivolumab_Riaz2017', 'PFS_Prog'],    
        
        ['Melanoma_Nivolumab_Riaz2017.On', 'OS_Naive'],
        ['Melanoma_Nivolumab_Riaz2017.On', 'PFS_Naive'],
        
        ['Melanoma_Nivolumab_Riaz2017.On', 'OS_Prog'],
        ['Melanoma_Nivolumab_Riaz2017.On', 'PFS_Prog'],    
        
        ['Melanoma_Nivolumab_Riaz2017.Diff', 'OS_Naive'],
        ['Melanoma_Nivolumab_Riaz2017.Diff', 'PFS_Naive'],
        
        ['Melanoma_Nivolumab_Riaz2017.Diff', 'OS_Prog'],
        ['Melanoma_Nivolumab_Riaz2017.Diff', 'PFS_Prog'],    
        
        ['Melanoma_PD1_Gide2019', 'OS'],
        ['Melanoma_PD1_Gide2019', 'PFS'],
        
        ['Melanoma_CTLA4+PD1_Gide2019', 'OS'],
        ['Melanoma_CTLA4+PD1_Gide2019', 'PFS'],
        
        ['Melanoma_PD1_Liu2019', 'OS_Naive'],
        ['Melanoma_PD1_Liu2019', 'PFS_Naive'],
        
        ['Melanoma_PD1_Liu2019', 'OS_Prog'],
        ['Melanoma_PD1_Liu2019', 'PFS_Prog'],    
        
        ['Melanoma_PD1_Cui2021', 'OS'],
        ['Melanoma_PD1_Cui2021', 'PFS'],
        
        ['CCRCC_Nivolumab_Braun2020', 'OS'],
        ['CCRCC_Nivolumab_Braun2020', 'PFS'],
        ['CCRCC_Nivolumab_Braun2020', 'irPFS'],
        
        ['CCRCC_Everolimus_Braun2020', 'OS'],
        ['CCRCC_Everolimus_Braun2020', 'PFS'],
        
        ['Urothelial_Atezolizumab_Mariathasan2018', 'OS'],
        
        ['NSCLC_Pembrolizumab_Rizvi2015', 'PFS'],
        
        ['mRCC_Atezolizumab_McDermott2018', 'PFS'],
        ['mRCC_Atezo+Bev_McDermott2018', 'PFS'],
        ['mRCC_Sunitinib_McDermott2018', 'PFS'],
        
        ['Melanoma_CTLA4_Snyder2014', 'OS'],
        ['Melanoma_CTLA4_Snyder2014.Post', 'OS'],
        
        ['Melanoma_ACT_Lauss2017', 'OS'],
        ['Melanoma_ACT_Lauss2017', 'PFS'],
        
        ['Hepatocellular_Atezolizumab_Finn2020', 'OS'],
        ['Hepatocellular_Atezolizumab_Finn2020', 'PFS'],
        
        ['Hepatocellular_Atezo+Bev_Finn2020', 'OS'],
        ['Hepatocellular_Atezo+Bev_Finn2020', 'PFS'],
        
        ['Hepatocellular_Sorafenib_Finn2020', 'OS'],
        ['Hepatocellular_Sorafenib_Finn2020', 'PFS'],
        
        ['GBM_PD1_Zhao2019', 'OS'],
        
        ['NSCLC_ICB_Ravi2023', 'OS'],
        ['NSCLC_ICB_Ravi2023', 'PFS'],
        
        ['NSCLC_ICB_Ravi2023.Adeno', 'OS'],
        ['NSCLC_ICB_Ravi2023.Adeno', 'PFS'],
        
        ['NSCLC_ICB_Ravi2023.Squamous', 'OS'],
        ['NSCLC_ICB_Ravi2023.Squamous', 'PFS'],
        
        ['PanCancer_ICB_Li2023', 'PFS'],
        ['PanCancer_ICB_Li2023', 'OS'],
        
        # filtered by low exon coverage
        #['Colorectal_ICB_Thibaudin2023', 'PFS'],
        #['Colorectal_ICB_Thibaudin2023', 'OS'],
        
        ['CCRCC_ICB_Miao2018', 'PFS'],
        ['CCRCC_ICB_Miao2018', 'OS'],
        
        ['Salivary_ICB_Vos2023', 'PFS'],
        
        ['SCLC_Atezolizumab_Nabet2024', 'OS'],
        ['SCLC_Atezolizumab_Nabet2024', 'PFS'],
        
        ['SCLC_Placebo_Nabet2024', 'OS'],
        ['SCLC_Placebo_Nabet2024', 'PFS'],
        
        ['Hepatocellular_ICB_Budhu2023', 'OS'],
        ['Hepatocellular_ICB_Budhu2023', 'PFS'],
        
        ['Esophageal_Durvalumab_Mamdani2021', 'RFS'],
        
        ['Urothelial_Atezolizumab_Snyder2017', 'PFS'],
        ['Urothelial_Atezolizumab_Snyder2017', 'OS'],
        
        ['Melanoma_ACT_Barras2024', 'PFS'],
        ['Melanoma_ACT_Barras2024', 'OS'],
        
        ['Melanoma_ACT_Barras2024.Post', 'PFS'],
        ['Melanoma_ACT_Barras2024.Post', 'OS'],
        
        ['RCC_Avelumab+Axitinib_Motzer2020', 'PFS'],
        ['RCC_Sunitinib_Motzer2020', 'PFS'],
        
        ['SCLC_Durvalumab_Roper2021', 'OS'],
        
        ['Urothelial_Atezo+Chemo_Hamidi2024-IMvigor130', 'OS'],
        ['Urothelial_Atezolizumab_Hamidi2024-IMvigor130','OS'],
        ['Urothelial_Chemotherapy_Hamidi2024-IMvigor130','OS'],
    
        ['Urothelial_Atezolizumab_Hamidi2024-IMvigor010','OS'],
        ['Urothelial_Surveillance_Hamidi2024-IMvigor010','OS'],
    
        ['Urothelial_Atezolizumab_Hamidi2024-IMvigor210','OS'],
    
        ['Urothelial_Atezolizumab_Hamidi2024-IMvigor211','OS'],
        ['Urothelial_Chemotherapy_Hamidi2024-IMvigor211','OS'],
    
        ['NSCLC_Atezolizumab_Patil2022-OAK','OS'],
        ['NSCLC_Atezolizumab_Patil2022-OAK','PFS'],
        
        ['NSCLC_Atezolizumab_Patil2022-POPLAR','OS'],
        ['NSCLC_Atezolizumab_Patil2022-POPLAR','PFS'],
        
        ['NSCLC_Docetaxel_Patil2022-OAK','OS'],
        ['NSCLC_Docetaxel_Patil2022-OAK','PFS'],
        
        ['NSCLC_Docetaxel_Patil2022-POPLAR','OS'],
        ['NSCLC_Docetaxel_Patil2022-POPLAR','PFS'],
        
        ['BiliaryTract_AtezoGemCis+Bev_Macarulla2024','OS'],
        ['BiliaryTract_AtezoGemCis+Bev_Macarulla2024','PFS'],
    
        ['BiliaryTract_AtezoGemCis+Plb_Macarulla2024','OS'],
        ['BiliaryTract_AtezoGemCis+Plb_Macarulla2024','PFS'],
        
        ['Urothelial_Atezolizumab_Powles2024','OS'],
        ['Urothelial_Atezolizumab_Powles2024','DFS'],
    
        ['Urothelial_Observation_Powles2024','OS'],
        ['Urothelial_Observation_Powles2024','DFS'],
        
        ['RCC_Atezo+Bev_Motzer2020-IMmotion151','PFS'],
        ['RCC_Sunitinib_Motzer2020-IMmotion151','PFS'],
    
        # partial data
        ['PanCancer_PD1_Prat2017', 'PFS'],
        
        ['HNSCC_ICB_Foy2022', 'PFS'],
        ['HNSCC_ICB_Foy2022', 'OS'],
        ['NSCLC_PD1_Foy2022', 'OS'],
    ]
    
    gz_postfix = re.compile('\.gz$')
    
    for title, survtype in IO_lst:
        fprefix = os.path.join(input_path, '*', title)
        
        data = fprefix + '.' + data_type + '.gz'
        
        # check exists
        data = glob(data)
        if len(data) == 0: continue
        
        assert len(data) == 1
        
        data = data.pop()
        fprefix = re.sub('\.' + data_type + '\.gz$', '', data)
        
        # also include the gene expression file, for potential CTL correction and T-cell dysfunction scores
        data_expression = fprefix + '.TPM.gz'
        
        if not os.path.exists(data_expression): data_expression = fprefix + '.expression.gz'
        
        # for some somatic data, expression will not exist
        if not os.path.exists(data_expression): data_expression = None
        
        clinical = os.path.join(os.path.dirname(fprefix), title.split('.')[0] + '.' + survtype)
        assert os.path.exists(clinical)
        
        # output is directly appended back to the data type input file
        output = gz_postfix.sub('', data) + '.response_' + survtype
        
        lst_survival.append([title + '.' + survtype, clinical, data, data_expression, output])
        


lst_response = []

def load_response_list(data_type):
    lst = []
    
    ################################################################
    # with binary outcome only
    fprefix = os.path.join(input_path, '*', 'Hepatocellular_PD1_Li2024')
    lst.append(['Hepatocellular_PD1_Li2024', fprefix, 'response'])
    
    fprefix = os.path.join(input_path, '*', 'Breast_Durvalumab+Olaparib_Pusztai2021')
    lst.append(['Breast_Durvalumab+Olaparib_Pusztai2021', fprefix, 'response'])
    
    fprefix = os.path.join(input_path, '*', 'Breast_Pembrolizumab_Wolf2022')
    #lst.append(['Breast_Pembrolizumab_Wolf2022', fprefix, 'response'])
    
    for subtype in ['Basal', 'Luminal']:
        lst.append(['Breast_Pembrolizumab_Wolf2022.' + subtype, fprefix, 'response_' + subtype])
    
    # filtered by low exon coverage
    #fprefix = os.path.join(input_path, '*', 'Esophageal_Atezolizumab_VanDenEnde2021')
    #lst.append(['Esophageal_Atezolizumab_VanDenEnde2021', fprefix, 'response'])
    
    fprefix = os.path.join(input_path, '*', 'HeadNeck_Pembrolizumab_Uppaluri2020')
    lst.append(['HeadNeck_Pembrolizumab_Uppaluri2020', fprefix, 'response'])
    
    fprefix = os.path.join(input_path, '*', 'NSCLC_Pembrolizumab_Lee2021')
    lst.append(['NSCLC_Pembrolizumab_Lee2021', fprefix, 'response'])
    
    fprefix = os.path.join(input_path, '*', 'Melanoma_MAGEA3_Montoya2013')
    lst.append(['Melanoma_MAGEA3_Montoya2013', fprefix, 'response'])
    
    for therapy in ['CTLA4', 'PD1']:
        title = 'Melanoma_%s_Roh2017' % therapy
        fprefix = os.path.join(input_path, '*', title)
        lst.append([title, fprefix, 'response'])
    
    for cancer in ['PanCancer', 'Melanoma', 'HeadNeck']:
        title = '%s_Pembrolizumab_Cristescu2018' % cancer
        fprefix = os.path.join(input_path, '*', title)
        lst.append([title, fprefix, 'response'])
    
    gz_postfix = re.compile('\.gz$')
    
    for title, fprefix, response_type in lst:
        data = fprefix + '.' + data_type + '.gz'
        
        data = glob(data)
        
        if len(data) == 0: continue
        
        if len(data) > 1:
            p = response_type.split('_').pop()
            
            lst = []
            
            for f in data:
                if f.find(p) >= 0:
                    lst.append(f)
            
            data = lst
        
        assert len(data) == 1
        data = data.pop()
        
        fprefix = re.sub('\.' + data_type + '\.gz$', '', data)
        
        clinical = os.path.join(os.path.dirname(fprefix), title.split('.')[0] + '.' + response_type)
        assert os.path.exists(clinical)
        
        data_expression = fprefix + '.TPM.gz'
        if not os.path.exists(data_expression): data_expression = fprefix + '.expression.gz'
        
        if not os.path.exists(data_expression): data_expression = None
        
        output = gz_postfix.sub('', data) + '.' + response_type
        
        lst_response.append([title, clinical, data, data_expression, output])



lst_RECIST = []

def load_RECIST_list(data_type):
    lst = []
    
    # as primary endpoints
    for title in [
        'NSCLC_PD1_Hu2023',
        'NSCLC_PD1_Yan2025',
        'mGC_Pembrolizumab_Kim2018',
        'PanCancer_Pembrolizumab_Yang2021',
        'Colorectal_ICB_Thibaudin2023',
        'Penile_Toripalimab_An2025',
        'Pediatric_Atezolizumab_Nabbi2023',
        ]:
        fprefix = os.path.join(input_path, '*', title)
        lst.append([title, fprefix, 'RECIST'])
    
    fprefix_template = os.path.join(input_path, '*', 'Melanoma_%s_Campbell2023')
    for treatment in ['CTLA4-to-PD1', 'PD1-to-CTLA4', 'CTLA4']:
        fprefix = fprefix_template % treatment
        lst.append([os.path.basename(fprefix), fprefix, 'RECIST'])
    
    gz_postfix = re.compile('\.gz$')
    
    for title, fprefix, response_type in lst:
        data = fprefix + '.' + data_type + '.gz'
        
        data = glob(data)
        if len(data) == 0: continue
        
        assert len(data) == 1
        data = data.pop()
        
        fprefix = re.sub('\.' + data_type + '\.gz$', '', data)
        
        clinical = os.path.join(os.path.dirname(fprefix), title.split('.')[0] + '.' + response_type)
        assert os.path.exists(clinical)
        
        data_expression = fprefix + '.TPM.gz'
        if not os.path.exists(data_expression):
            data_expression = fprefix + '.expression.gz'
        
        if not os.path.exists(data_expression): data_expression = None
        
        output = gz_postfix.sub('', data) + '.response_' + response_type
        lst_RECIST.append([title, clinical, data, data_expression, output])




def load_clinical_data_pairs(f, data, data_expression = None, purity_covariate = None):
    clinical = pandas.read_csv(f, sep='\t', index_col=0)
    clinical.index = clinical.index.astype(str)
    
    flag_ICGC = os.path.dirname(data).find('ICGC') > 0
    
    data = pandas.read_csv(data, sep='\t', index_col=0)
    
    flag_TCGA = (data.columns[0].find('TCGA-') == 0)
    flag_TARGET = (data.columns[0].find('TARGET-') == 0)
    
    # if TCGA names
    if flag_TCGA or flag_TARGET:
        # only look at non-normal samples
        data = data.loc[:, [v.split('.').pop() != 'Normal' for v in data.columns]]
        
        # merge by patient identities
        flag = ['-'.join(v.split('-')[:3]) for v in data.columns]
    
    elif flag_ICGC:
        # merge by patient identities
        flag = [v.split('.')[0] for v in data.columns]
    
    else:
        flag = None
    
    if flag is not None:
        data.columns = flag
        data = data.T.groupby(data.columns).median().T
        data = data.loc[(data == 0).mean(axis=1) < 1]
    
    common = clinical.index.intersection(data.columns)
    
    flag_purity_covariate = False
    
    if purity_covariate is not None:
        # correction pair: data, field name
        if 'purity' in [v.lower() for v in clinical.columns]:
            print('purity already exist in clinical covariates')
        
        else:
            purity_f = '.'.join(f.split('.')[:-1]) + '.Purity.gz'
            
            if os.path.exists(purity_f):
                print('use original purity values')
                purity_covariate = pandas.read_csv(purity_f, sep='\t', index_col=0)['Purity']
                purity_covariate.index = purity_covariate.index.astype(str)
            
            else:
                purity_covariate, title = purity_covariate
                
                print('use computed purity values', purity_covariate)
                purity_covariate = pandas.read_csv(purity_covariate, sep='\t', index_col=0)
                purity_covariate = purity_covariate.loc[title]
                
            common = common.intersection(purity_covariate.index)
            
            flag_purity_covariate = True    
    
    if data_expression is None:
        clinical, data = clinical.loc[common], data.loc[:, common]
        clinical = clinical.loc[:, clinical.abs().sum() > 0]
        
        if flag_purity_covariate:
            clinical[purity_covariate.name] = purity_covariate.loc[common]
        
        return clinical, data
    
    else:
        data_expression = pandas.read_csv(data_expression, sep='\t', index_col=0)
        
        if flag_TCGA or flag_TARGET:
            flag = ['-'.join(v.split('-')[:3]) for v in data_expression.columns]
        
        elif flag_ICGC:
            # merge by patient identities
            flag = [v.split('.')[0] for v in data_expression.columns]
        
        else:
            flag = None
        
        if flag is not None:
            data_expression = data_expression.T.groupby(flag).median().T
            data_expression = data_expression.loc[(data_expression == 0).mean(axis=1) < 1]
        
        common = common.intersection(data_expression.columns)
        
        clinical, data, data_expression = clinical.loc[common], data.loc[:, common], data_expression.loc[:, common]
        clinical = clinical.loc[:, clinical.abs().sum() > 0]
        
        if flag_purity_covariate:
            clinical[purity_covariate.name] = purity_covariate.loc[common]
        
        return clinical, data, data_expression





def compute_survival_associations(flag_CTL_correction, best_threshold_margin, cnt_thres = 10, purity_covariate = None):
    # purity_covariate: postfix of background covariate estimated value
    inx, Nnode = int(sys.argv[1]), int(sys.argv[2])
    
    if len(lst_survival) != Nnode:
        sys.stderr.write('Please use ' + str(len(lst_survival)) + ' as N.\n')
        sys.exit(1)
    
    title, clinical, data_input, data_expression, output = lst_survival[inx]
    print('process', title)
    
    if purity_covariate is not None:
        output += '.purity_corrected'
        
        f_template = '.'.join(data_input.replace('.gz', '').split('.')[:-1]) + '.%s.' + purity_covariate[0]
        
        fpurity = glob(f_template % 'TPM')
        fpurity.extend(glob(f_template % 'expression'))
        
        if len(fpurity) == 0:
            print('Cannot find purity covariate %s for %s' % (purity_covariate[0], data_input))
            return
        
        assert len(fpurity) == 1
        purity_covariate = [fpurity[0], purity_covariate[1]]
    
    cmd_prefix = os.path.join(base_path, 'src', 'survival_associations.R')
    
    # always load all data, as there could be many modifications
    if flag_CTL_correction:
        # no expression file, just return
        if data_expression is None: return
        clinical, data, data_expression = load_clinical_data_pairs(clinical, data_input, data_expression, purity_covariate = purity_covariate)
    
        output += '.CTL_corrected'
        
        CTL_set = ['CD8A', 'CD8B', 'GZMA', 'GZMB', 'PRF1']
        
        if len(data_expression.index.intersection(CTL_set)) <= len(CTL_set)/2:
            sys.stderr.write('Insufficient CD8 T genes for %s\n' % title)
            return
        
        clinical['CTL'] = data_expression.reindex(CTL_set).dropna().median()
    
    else:
        # no need to load the gene expression part, which may undermine the sample number
        clinical, data = load_clinical_data_pairs(clinical, data_input, purity_covariate = purity_covariate)
    
    if clinical.shape[0] < cnt_thres:
        print('skip %s by low sample size %d.' % (os.path.basename(output), clinical.shape[0]))
        return
    
    #if os.path.exists(output): return
    
    # some fix here for special dataset
    #if title == 'ORCA-IN.OS':
    #    # after overlapping with expression, this will cause too few variations
    #    clinical.drop('Stage', axis=1, inplace=True)
    
    # if is Mutation file, correct for base mutation burden
    if output.find('.Mutation.response') > 0:
        # log2 transform mutation burden
        clinical['mutation burden'] = numpy.log2(data.loc['mutation_burden'] + 1)    
        data.drop('mutation_burden', inplace=True)
    
    if output.find('.CNA.response') > 0:
        clinical['Aneuploidy'] = data.loc['Aneuploidy']    
        data.drop('Aneuploidy', inplace=True)
    
    # MiXCR may have NA
    data.fillna(0, inplace=True)
    
    # as overlap operations happend before, check the data matrix again
    data = data.loc[(data == 0).mean(axis=1) < 1]
    
    if data.shape[0] == 0:
        print('skip %s by all zero values.' % os.path.basename(output))
        return
    
    # write download modified clinical, data
    data.to_csv(output + '.data_temp', sep='\t', index_label=False)    
    clinical.to_csv(output + '.clinical', sep='\t', index_label=False)
    
    cmd = [cmd_prefix, output + '.data_temp', output + '.clinical', output]
        
    if best_threshold_margin is not None: cmd.append(str(best_threshold_margin))
    
    os.system(' '.join(cmd))
    
    if clinical.shape[1] > 2:
        # has some factors, generate the clinical along survival effect
        md5_prev = hashlib.md5(open(output + '.clinical','rb').read()).hexdigest()
        
        cmd = [os.path.join(base_path, 'src', 'coxph.R'), output + '.clinical', output + '.clinical']
        os.system(' '.join(cmd))
        
        md5 = hashlib.md5(open(output + '.clinical','rb').read()).hexdigest()
        if md5 == md5_prev:
            print('delete clinical file for %s' % output)
            os.remove(output + '.clinical')
        
    else:
        os.remove(output + '.clinical')
    
    for postfix in ['data_temp']: os.remove(output + '.' + postfix)




def compute_binary_response_associations(flag_CTL_correction, purity_covariate=None):
    # more complicated version, by logistic regression
    
    for title, response, data_input, data_expression, output in lst_response:
        
        if purity_covariate is not None:
            output += '.purity_corrected'
            
            f_template = '.'.join(data_input.replace('.gz', '').split('.')[:-1]) + '.%s.' + purity_covariate[0]
            
            fpurity = glob(f_template % 'TPM')
            fpurity.extend(glob(f_template % 'expression'))
            
            if len(fpurity) == 0:
                print('Cannot find purity covariate %s for %s' % (purity_covariate[0], data_input))
                return
            
            assert len(fpurity) == 1
            purity_covariate = [fpurity[0], purity_covariate[1]]
        
        if flag_CTL_correction:
            output += '.CTL_corrected'
                
            # no CTL expression for correction
            if data_expression is None:
                print('Cannot find expression data for %s' % data_input)
                return
            
            response, data, data_expression = load_clinical_data_pairs(response, data_input, data_expression, purity_covariate=purity_covariate)
                            
            if title.find('Mouse') >= 0:
                CTL_set = ['CD8a', 'Cd8b1', 'Gzma', 'Gzmb', 'Prf1']
            else:
                CTL_set = ['CD8A', 'CD8B', 'GZMA', 'GZMB', 'PRF1']
                
            CTL = data_expression.reindex(CTL_set).dropna().median()
                
            response['CTL'] = CTL
            
        else:
            response, data = load_clinical_data_pairs(response, data_input, purity_covariate=purity_covariate)
            
            
        if output.find('.Mutation.response') > 0:
            response['mutation burden'] = numpy.log2(data.loc['mutation_burden'] + 1)
            data.drop('mutation_burden', inplace=True)
            
        if output.find('.CNA.response') > 0:
            response['Aneuploidy'] = data.loc['Aneuploidy']
            data.drop('Aneuploidy', inplace=True)
            
        # has clinical cov
        if response.shape[1] > 1:
            response['Const'] = 1
                
            y = CytoSig.dataframe_to_array(response.iloc[:, 0].astype(float))
            X = CytoSig.dataframe_to_array(response.iloc[:, 1:].astype(float))
        
            try:
                _, _, z, p = data_significance.logit(X, y, 1e-5, 1, 1000, 1, 0)
            except ArithmeticError:
                sys.stderr.write('Clinical file only is even not regressable. Exit\n')
                return
                
            result = pandas.DataFrame([z, p], index=['z', 'p'], columns=response.columns[1:]).transpose()
            result.to_csv(output + '.clinical', sep='\t', index_label=False)
                
            response.drop('Const', axis=1, inplace=True)
        
        
        print('process', title, response.shape[0], 'samples')
            
        # fill NA here to allow the computation
        data.fillna(0, inplace=True)
        
        if response.shape[1] == 1:
            print('only response is involved, use ranksum')
                
            response = response.iloc[:, 0].astype(bool)
                
            stat = data.apply(
                lambda arr: pandas.Series(
                    stats.ranksums(arr.loc[response], arr.loc[~response]), index=['z', 'p'], name=arr.name)
                , axis=1)
                
        else:
            print('multivariate, use logistic regression', response.iloc[:, 0].sum(), response.iloc[:, 0].mean())
        
            response.iloc[:, 0].to_csv(output + '.y', sep='\t', index_label=False)
            response.iloc[:, 1:].to_csv(output + '.B', sep='\t', index_label=False)
            data.transpose().to_csv(output + '.X', sep='\t', index_label=False)
        
            cmd = os.path.join('logis_batch -cntthres 3 -verbose 0')
            os.system(' '.join([cmd, '-B', output + '.B', '-Y', output + '.y', '-X', output + '.X', '-out', output]))
        
            z = pandas.read_csv(output + '.zscore', sep='\t', index_col=0).iloc[:, 0]
            z.name = 'z'
        
            p = pandas.read_csv(output + '.pvalue', sep='\t', index_col=0).iloc[:, 0]
            p.name = 'p'
        
            stat = pandas.concat([z, p], axis=1, join='inner')
            for postfix in ['zscore', 'pvalue', 'coef', 'B', 'X', 'y']: os.remove(output + '.' + postfix)
                
            response = response.iloc[:, 0].astype(bool)
            
            
        stat['mean.value'] = data.loc[stat.index].mean(axis=1)
        stat['N'] = response.shape[0]
        stat['FDR'] = multipletests(stat['p'], method='fdr_bh')[1]
            
        stat.to_csv(output, sep='\t', index_label=False)
        
        if not flag_CTL_correction:
            # binary outcome, add AUC
            merge = []
            for title, arr in data.iterrows():
                if title in ['TIDE', 'TIDE.exclusion', 'CytoSig.TGFB1', 'Signature.TGFb', 'Aneuploidy']:
                    arr = -arr
                
                fpr, tpr, _ = metrics.roc_curve(response, arr)
                AUC = metrics.auc(fpr, tpr)
                p = stats.mannwhitneyu(arr.loc[response], arr.loc[~response])[1]
                
                merge.append(pandas.Series([AUC, p], index=['AUC', 'p'], name=title))
                #lst.append([title, metrics.auc(fpr, tpr), fpr, tpr])
            
            stat = pandas.concat(merge, axis=1, join='inner').T
            stat['mean.value'] = data.loc[stat.index].mean(axis=1)
            stat['N'] = response.shape[0]
            stat['FDR'] = multipletests(stat['p'], method='fdr_bh')[1]
                
            stat.to_csv(output + '_AUC', sep='\t', index_label=False)



def compute_RECIST_response_associations(flag_CTL_correction, cnt_thres = 10, purity_covariate=None):
    warnings.filterwarnings("error")
    
    for title, clinical, data_input, data_expression, output in lst_RECIST:
        print(title)
        
        if purity_covariate is not None:
            output += '.purity_corrected'
            
            f_template = '.'.join(data_input.replace('.gz', '').split('.')[:-1]) + '.%s.' + purity_covariate[0]
        
            fpurity = glob(f_template % 'TPM')
            fpurity.extend(glob(f_template % 'expression'))
        
            if len(fpurity) == 0:
                print('Cannot find purity covariate %s for %s' % (purity_covariate[0], data_input))
                return
            
            assert len(fpurity) == 1
            purity_covariate_i = [fpurity[0], purity_covariate[1]]
        else:
            purity_covariate_i = None
    
        if flag_CTL_correction:
            # expression not existing, just jump
            if data_expression is None: continue
            clinical, data, data_expression = load_clinical_data_pairs(clinical, data_input, data_expression, purity_covariate=purity_covariate_i)
            output += '.CTL_corrected'
            
            CTL = data_expression.reindex(['CD8A', 'CD8B', 'GZMA', 'GZMB', 'PRF1']).dropna().median()
            clinical['CTL'] = CTL
        
        else:
            clinical, data = load_clinical_data_pairs(clinical, data_input, purity_covariate=purity_covariate_i)
        
        # jump too few samples
        if clinical.shape[0] < cnt_thres: continue
        
        # if is mutation
        if output.find('.Mutation.response') > 0:
            clinical['mutation burden'] = numpy.log2(data.loc['mutation_burden'] + 1)
            data.drop('mutation_burden', inplace=True)
        
        if output.find('.CNA.response') > 0:
            clinical['Aneuploidy'] = data.loc['Aneuploidy']    
            data.drop('Aneuploidy', inplace=True)

        clinical['Const'] = 1.0
        
        data.fillna(0, inplace=True)
        
        # analyze base clinical factors first
        result = sm.OLS(clinical.iloc[:, 0], clinical.iloc[:, 1:]).fit()
        result = pandas.concat([result.tvalues, result.pvalues], axis=1)
        result.columns = ['t', 'p']
        result.to_csv(output + '.clinical', sep='\t', index_label=False)
        
        clinical['Target'] = 1.0
        
        print(clinical.shape[0], 'samples')
    
        y = clinical.iloc[:, 0]
        X = clinical.iloc[:, 1:]
        
        merge = []
        
        for pivot, arr in data.iterrows():
            #if pivot.split('@')[-1] not in ['FIBP', 'FCMR']: continue
            if arr.std() == 0: continue
            
            X.loc[:, 'Target'] = arr.astype(float)
            
            try:
                result = sm.OLS(y, X).fit()
                
                merge.append(
                    pandas.Series([result.tvalues['Target'], result.pvalues['Target']],
                        index=['t', 'p'], name=pivot)
                    )
            
            except (RuntimeWarning, MissingDataError):
                continue
        
        stat = pandas.concat(merge, axis=1).transpose()
        
        stat['mean.value'] = data.loc[stat.index].mean(axis=1)
        stat['N'] = clinical.shape[0]
        stat['FDR'] = multipletests(stat['p'], method='fdr_bh')[1]
        
        stat.to_csv(output, sep='\t', index_label=False)
        print(stat.shape)




def main():
    # don't correct for CTL infiltration effects for all analysis
    flag_CTL_correction = False
    
    # best separation threshold selection margin
    best_threshold_margin = 5
    
    #flag_CTL_correction, best_threshold_margin = True, None
    
    # load data files
    for data_type in [
        # gene expression
        'expression', 'TPM',
        # isoform expression mapped from raw data
        'isoform_TPM',
        # somatic alterations
        'Mutation' , 'CNA',
        # DNA promoter methylation
        'TSS200_methylation',
        ]:
        load_survival_list(data_type)
        load_response_list(data_type)
        load_RECIST_list(data_type)    
    
    if len(sys.argv) == 3:
        compute_survival_associations(flag_CTL_correction, best_threshold_margin)
    else:
        compute_binary_response_associations(flag_CTL_correction)
        compute_RECIST_response_associations(flag_CTL_correction)
    
    return 0

if __name__ == '__main__': main()
