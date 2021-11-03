import numpy as np
import os
import argparse
from dohmmscan import hmmscan
from gb2faa import gb2faa
from viterbi import viterbi
from sklearn.metrics import roc_curve, auc
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

B= [[], []]
A = []
DICT = {}
start_p = []

def mkdir(path):
    mkdircmd = 'mkdir -p %s' % path
    os.system(mkdircmd)

def init_model(acc_or_pred = 1):
    with open('%s/emission_prob'%model_path, 'r') as femiss:
        for line in femiss.readlines():
            B[0].append(float(line.split('\t')[1]))
            B[1].append(float(line.split('\t')[2]))
    if acc_or_pred:
        B[0].append(1.)
        B[1].append(1.)
    else:
        B[0].append(0.)
        B[1].append(0.)
    with open('%s/pfam_DICT'%model_path, 'r') as fdict:
        for line in fdict.readlines():
            DICT[line.split('\t')[0]] = int(line.split('\t')[1])
    with open('%s/transition_prob'%model_path, 'r') as ftrans:
        for line in ftrans.readlines():
            temp = []
            for i in line.strip().split('\t'):
                temp.append(float(i))
            A.append(temp)
    with open('%s/start_prob'%model_path, 'r') as fstart:
        for line in fstart.readlines():
            for i in line.strip().split('\t'):
                start_p.append(float(i))
    mkdir(test_dir)
    mkdir(test_faa_path)
    mkdir(test_hmmscan_path)
    mkdir(result_path)
    mkdir(each_phage_result_path)

def testsetall(testphage_list, evalue):
    score = []
    true_list = []
    fout = open('%s/all_test_res_table.txt'%result_path, 'w')
    fout.write('PhageId' + '\t' + 'ProteinNum' + '\t' + 'ActTailNum' + '\t' + 'PredictTailNum' + '\t' + 'TPR' + '\t' + 'FNR' + '\t' + 'TNR' + '\t' + 'FPR' + '\t' + 'Observation' + '\t' + 'Predict' + '\n')

    ftest = open(testphage_list, 'r')
    for line in ftest.readlines():

        observation = []
        phageid = line.split('\t')[0].split('.')[0].strip()
        phageid_withver = line.split('\t')[0].strip()
        print('----------%s is now testing----------'%phageid)

        phagepath = os.path.join(test_faa_path, phageid_withver) + '.faa'
        gb2faa(phageid_withver, test_faa_path)
        hmmscan(phagepath, test_hmmscan_path, hmmdb, evalue)

        prot_content_list = []
        prot_list = []
        prot_start_list = []
        prot_end_list = []
        prot_def_list = []
        prot_size_list = []
        count_list = []
        count = 0
        with open(phagepath, 'r') as fphage:
            for line in fphage.readlines():
                if '>' in line:
                    count += 1
                    prot_content = line.split('>ref|')[-1].strip()
                    protid = line.split('\t')[0].split('|')[1]
                    prot_start = line.split('\t')[0].split('|')[3].split(':')[0].split('[')[-1]
                    prot_end = line.split('\t')[0].split('|')[3].split(':')[1].split(']')[0]
                    prot_def = line.split('\t')[0].split('|')[4]
                    prot_size = (int(prot_end) - int(prot_start)) // 3
                    prot_content_list.append(prot_content)
                    prot_list.append(protid)
                    prot_start_list.append(prot_start)
                    prot_end_list.append(prot_end)
                    prot_def_list.append(prot_def)
                    prot_size_list.append(prot_size)
                    count_list.append(str(count))
                    observation.append(-1)

        with open('%s/%snontail.cpt'%(test_hmmscan_path, phageid), 'r') as fnontail:
            nontailcontext = fnontail.readlines()
            if nontailcontext:
                for line in nontailcontext:
                    thisprotid = line.split('\t')[0].split('|')[1]
                    count = prot_list.index(thisprotid)
                    pfamid = line.split('\t')[1]
                    if pfamid in DICT.keys():
                        observation[count] = DICT[pfamid]
        with open('%s/%stail.cpt'%(test_hmmscan_path, phageid), 'r') as ftail:
            tailcontext = ftail.readlines()
            if tailcontext:
                for line in tailcontext:
                    thisprotid = line.split('\t')[0].split('|')[1]
                    count = prot_list.index(thisprotid)
                    pfamid = line.split('\t')[1]
                    if pfamid in DICT.keys():
                        observation[count] = DICT[pfamid]
            else:
                print('%s no tail in def' % phageid)

        res, F = viterbi(A, B, start_p, observation)
        for i in range(len(res)):
            true_list.append(int(res[i]))
            each_score = F[0][i] - F[1][i]
            score.append(each_score)
        protcount = 0
        act_tailcount = 0
        act_nontailcount = 0
        pred_tailcount = 0
        pred_nontailcount = 0
        TPcount, FPcount, TNcount, FNcount = 0, 0, 0, 0
        for i in range(len(observation)):
            protcount += 1
            if 0 <= observation[i] <= 839:
                act_tailcount += 1
                tail_act = True
            else:
                act_nontailcount += 1
                tail_act = False
            if res[i] == 1:
                pred_tailcount += 1
                tail_pred = True
            else:
                pred_nontailcount += 1
                tail_pred = False

            if tail_act and tail_pred:
                TPcount += 1
            elif not tail_act and tail_pred:
                FPcount += 1
            elif not tail_act and not tail_pred:
                TNcount += 1
            else:
                FNcount += 1
    
        TPR, FNR, TNR, FPR = 0, 0, 0, 0
        try:
            TPR = TPcount / act_tailcount
        except:
            TPR = 0
        try:
            FPR = FPcount / act_nontailcount
        except:
            FPR = 0
        try:
            FNR = FNcount / act_tailcount
        except:
            FNR = 0
        try:
            TNR = TNcount / act_nontailcount
        except:
            TNR = 0

        observation = list(map(str, observation))
        res = list(map(str, res))

        fout_phage = open('%s/%s_prot_result_table.txt'%(each_phage_result_path, phageid), 'w')
        fout_phage.write('PhageContent' + '\t' + 'PhageId' + '\t' + 'ProteinId' + '\t' + 'ProteinSize' + '\t' + 'ProteinCount' + '\t' + 'ProteinStartIndex' + '\t' + 'ProteinEndIndex' + '\t' + 'ProteinDef' + '\t' + 'TailOrNotTail' + '\n' )
        for i in range(count):
            fout_phage.write(prot_content_list[i] + '\t' + phageid_withver + '\t' + prot_list[i] + '\t' + str(prot_size_list[i]) + '\t' + str(i) + '\t' + str(prot_start_list[i]) + '\t' + str(prot_end_list[i])  + '\t' + prot_def_list[i] + '\t' + res[i] + '\n')
        fout_phage.close()

        obs = ','.join(observation)
        pred = ','.join(res)
        fout.write(phageid_withver + '\t' + str(protcount) + '\t' + str(act_tailcount) + '\t' + str(pred_tailcount) + '\t' + str(TPR) + '\t' + str(FNR) + '\t' + str(TNR) + '\t' + str(FPR) + '\t' + obs + '\t' + pred + '\n')

    ftest.close()
    fout.close()

    fpr, tpr, thresholds = roc_curve(true_list, score)
    AUC = auc(fpr, tpr)
    plt.plot(fpr,tpr,label='HMM model AUC %0.2f' % AUC, color='blue', lw = 2)
    plt.plot([0, 1], [0, 1], color='black', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating Curve')
    plt.legend(loc="lower right")
    plt.savefig('%s/testset_ROC.jpg'%result_path)

def testone(inputphage, evalue):
    if not os.path.exists(inputphage):
        if '.' in inputphage:
            phageid = inputphage.split('.')[0]
        else:
            phageid = inputphage
        phage_path = os.path.join(test_faa_path, phageid) + '.faa'
        gb2faa(phageid, test_faa_path)
    else:
        phage_path = inputphage
        phageid = inputphage.split('/')[-1].split('.')[0]

    observation = []
    hmmscan(phage_path, test_hmmscan_path, hmmdb, evalue)

    prot_content_list = []
    prot_list = []
    prot_start_list = []
    prot_end_list = []
    prot_def_list = []
    prot_size_list = []
    count_list = []
    count = 0
    with open(phage_path, 'r') as fphage:
        for line in fphage.readlines():
            if '>' in line:
                count += 1
                prot_content = line.split('>ref|')[-1].strip()
                protid = line.split('\t')[0].split('|')[1]
                prot_start = line.split('\t')[0].split('|')[3].split(':')[0].split('[')[-1]
                prot_end = line.split('\t')[0].split('|')[3].split(':')[1].split(']')[0]
                prot_def = line.split('\t')[0].split('|')[4]
                prot_size = (int(prot_end) - int(prot_start)) // 3
                prot_content_list.append(prot_content)
                prot_list.append(protid)
                prot_start_list.append(prot_start)
                prot_end_list.append(prot_end)
                prot_def_list.append(prot_def)
                prot_size_list.append(prot_size)
                count_list.append(str(count))
                observation.append(-1)

    with open('%s/%snontail.cpt'%(test_hmmscan_path, phageid), 'r') as fnontail:
        nontailcontext = fnontail.readlines()
        if nontailcontext:
            for line in nontailcontext:
                thisprotid = line.split('\t')[0].split('|')[1]
                count = prot_list.index(thisprotid)
                pfamid = line.split('\t')[1]
                if pfamid in DICT.keys():
                    observation[count] = DICT[pfamid]
    with open('%s/%stail.cpt'%(test_hmmscan_path, phageid), 'r') as ftail:
        tailcontext = ftail.readlines()
        if tailcontext:
            for line in tailcontext:
                thisprotid = line.split('\t')[0].split('|')[1]
                count = prot_list.index(thisprotid)
                pfamid = line.split('\t')[1]
                if pfamid in DICT.keys():
                    observation[count] = DICT[pfamid]
        else:
            print('%s no tail in def' % phageid)

    res, F = viterbi(A, B, start_p, observation)
    protcount = 0
    act_tailcount = 0
    act_nontailcount = 0
    pred_tailcount = 0
    pred_nontailcount = 0
    TPcount, FPcount, TNcount, FNcount = 0, 0, 0, 0
    for i in range(len(observation)):
        protcount += 1
        if 0 <= observation[i] <= 839:
            act_tailcount += 1
            tail_act = True
        else:
            act_nontailcount += 1
            tail_act = False
        if res[i] == 1:
            pred_tailcount += 1
            tail_pred = True
        else:
            pred_nontailcount += 1
            tail_pred = False

        if tail_act and tail_pred:
            TPcount += 1
        elif not tail_act and tail_pred:
            FPcount += 1
        elif not tail_act and not tail_pred:
            TNcount += 1
        else:
            FNcount += 1
    
    TPR, FNR, TNR, FPR = 0, 0, 0, 0
    try:
        TPR = TPcount / act_tailcount
    except:
        TPR = 0
    try:
        FPR = FPcount / act_nontailcount
    except:
        FPR = 0
    try:
        FNR = FNcount / act_tailcount
    except:
        FNR = 0
    try:
        TNR = TNcount / act_nontailcount
    except:
        TNR = 0

    observation = list(map(str, observation))
    res = list(map(str, res))

    fout_phage = open('%s/%s_prot_result_table.txt'%(result_path, phageid), 'w')
    fout_phage.write('PhageContent' + '\t' + 'PhageId' + '\t' + 'ProteinId' + '\t' + 'ProteinSize' + '\t' + 'ProteinCount' + '\t' + 'ProteinStartIndex' + '\t' + 'ProteinEndIndex' + '\t' + 'ProteinDef' + '\t' + 'TailOrNotTail' + '\n' )
    for i in range(count):
        fout_phage.write(prot_content_list[i] + '\t' + phageid + '\t' + prot_list[i] + '\t' + str(prot_size_list[i]) + '\t' + str(i) + '\t' + str(prot_start_list[i]) + '\t' + str(prot_end_list[i])  + '\t' + prot_def_list[i] + '\t' + res[i] + '\n')
    fout_phage.close()

    fout = open('%s/%s_result_table.txt'%(each_phage_result_path, phageid), 'w')
    fout.write('PhageId' + '\t' + 'ProteinNum' + '\t' + 'ActTailNum' + '\t' + 'PredictTailNum' + '\t' + 'TPR' + '\t' + 'FNR' + '\t' + 'TNR' + '\t' + 'FTR' + '\t' + 'Observation' + '\t' + 'Predict' + '\n')
    obs = ','.join(observation)
    pred = ','.join(res)
    fout.write(phageid + '\t' + str(protcount) + '\t' + str(act_tailcount) + '\t' + str(pred_tailcount) + '\t' + str(TPR) + '\t' + str(FNR) + '\t' + str(TNR) + '\t' + str(FPR) + '\t' + obs + '\t' + pred + '\n')
    fout.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A HMM Model For Tail ')
    parser.add_argument('-i', '--input', required='True', help='abs_path or input phage faa or phage list, \ninput support phageid or phage faa file')
    parser.add_argument('-o', '--output', required='True', help='output abs_path')
    parser.add_argument('--accurate--mode', action='store_true', default=False, help='Change predict mode to accurate')
    parser.add_argument('--hmmscan--evalue', type=float, default=10.0, help='Specifies a evalue threshold for the hmmscan(default 10.0)')
    parser.add_argument('--phagelist--mode', action='store_true', default=False, help='Specifies input a phage faa or phage list to test, default input a pahge faa')
    args = parser.parse_args()

    result_path = args.output

    each_phage_result_path = '%s/each_phage_result'%result_path
    root_path = os.path.abspath("..")
    model_path = '%s/hmmmodel'%root_path
    train_dir = '%s/train_data'%root_path
    train_phage_list = '%s/train_phage_list.txt'%train_dir
    test_phage_list = '%s/test_phage_list.txt'%root_path
    test_dir = '%s/test_data'%result_path
    test_faa_path = '%s/test_phage_faa'%test_dir
    test_hmmscan_path = '%s/test_hmmscan_res'%test_dir
    hmmdb = '%s/dbs'%root_path

    if args.accurate__mode:
        acc_or_pred = 0
        init_model(acc_or_pred)
    else:
        init_model()

    if args.phagelist__mode:
        test_phage_list = args.input
        testsetall(test_phage_list, args.hmmscan__evalue)
    else:
        phagefaa = args.input
        testone(phagefaa, args.hmmscan__evalue)