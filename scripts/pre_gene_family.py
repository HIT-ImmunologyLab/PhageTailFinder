import os

def mkdir(path):
    mkdircmd = 'mkdir -p %s' % path
    os.system(mkdircmd)

genefamily_path = '/home/yangh/PTF/result/gene_family_result'
test_phage_path = '/home/yangh/PTF/result/test_data/test_phage_faa'
mkdir(genefamily_path)

fout_table = open('%s/test_tail_table.txt'%genefamily_path, 'w')
fout_faa = open('%s/test_tail.faa'%genefamily_path, 'w')

phage_family_dict = {}
for i in os.listdir('/home/yangh/PhageTailFinder/all_tailedphage_data/tailed_phage'):
    family_path = os.path.join('/home/yangh/PhageTailFinder/all_tailedphage_data/tailed_phage', i)
    for j in os.listdir(family_path):
        phage_family_dict[j.split('.faa')[0]] = i

fout_table.write('protein_id' + '\t' + 'protein_size' + '\t' + 'phage_id' + '\t' + 'protein_def' + '\t' + 'protein_start_end' + '\t' + 'phage_family' + '\t' + 'count' + '\n')

for i in os.listdir('/home/yangh/PTF/result/each_phage_result'):
    each_phage_table_path = os.path.join('/home/yangh/PTF/result/each_phage_result', i)
    with open(each_phage_table_path, 'r') as fin:
        count = 0
        for line in fin.readlines()[1:]:
            prot_content = line.split('\t')[0]
            phage_id = line.split('\t')[1]
            prot_id = line.split('\t')[2]
            prot_size = line.split('\t')[3]
            prot_count = line.split('\t')[4]
            prot_start = line.split('\t')[5]
            prot_end = line.split('\t')[6]
            prot_def = line.split('\t')[7]
            tail_flag = line.split('\t')[8].strip()
            phage_family = phage_family_dict[phage_id]
            with open('%s/%s.faa'%(test_phage_path, phage_id), 'r') as fphage:
                phage_content = fphage.readlines()
            if tail_flag == '1':
                count += 1
                fout_faa.write('>' + phage_content[int(prot_count)*2].split('>ref|')[-1])
                fout_faa.write(phage_content[int(prot_count)*2+1])
                fout_table.write(prot_content + '\t' + prot_size + '\t' + phage_id + '\t' + prot_def + '\t' + prot_start + ':' + prot_end + '\t' + phage_family + '\t' + str(count) + '\n')
fout_table.close()
fout_faa.close()