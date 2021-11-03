import os

def getprotinfo(inputfile, hmmscan_res):
    faa_id = inputfile.split('/')[-1].split('.')[0]
    f1 = open(hmmscan_res + '/' + faa_id + '.xml', 'r')
    genequery = {}
    for line in f1:
        if line[0] == '#': continue
        linetab = line.split()
        if linetab[2] not in genequery.keys():
            genequery[linetab[2]] = []
            genequery[linetab[2]].append(linetab[0])  # 0 tail protein name
            genequery[linetab[2]].append(linetab[2])  # 1 query name
            genequery[linetab[2]].append(linetab[4])  # 2 e-value
            genequery[linetab[2]].append(linetab[5])  # 3 score
            genequery[linetab[2]].append(linetab[6])  # 4 bias
            continue
    f2 = open(hmmscan_res + '/' + faa_id + '.cpt', 'w')
    genequery = sorted(genequery.items())

    for item in genequery:
        # print(item)
        f2.write('%s\t%s\t%s\t%s\t%s\n' % (item[0], item[1][0], item[1][2], item[1][3], item[1][4]))
    f1.close()
    f2.close()
    print(faa_id + '  xml2cpt success')

def hmmscan(inputfile, hmmscan_res, hmmdb, evalue = 10):
    inputid = inputfile.split('/')[-1].split('.')[0]
    if os.path.exists('%s/%stail.cpt'%(hmmscan_res, inputid)):
        return 0
    print('%s now doing hmmscan'%inputid)
    tailhmmprofile = hmmdb + '/tail_pfam'
    nontailhmmprofile = hmmdb + '/nontail_pfam'
    tailhmmscancmd = 'hmmscan --tblout %s/%s.xml %s %s > /dev/null' % (hmmscan_res, inputid+'tail', tailhmmprofile, inputfile)
    nontailhmmscancmd = 'hmmscan -E %s --tblout %s/%s.xml %s %s > /dev/null' % (evalue, hmmscan_res, inputid+'nontail', nontailhmmprofile, inputfile)
    os.system(tailhmmscancmd)
    os.system(nontailhmmscancmd)
    getprotinfo(inputid+'tail', hmmscan_res)
    getprotinfo(inputid+'nontail', hmmscan_res)