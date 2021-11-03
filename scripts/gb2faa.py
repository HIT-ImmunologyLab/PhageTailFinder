import os
from Bio import SeqIO
from Bio import Entrez
import ssl
Entrez.email = '979268136@qq.com'
ssl._create_default_https_context = ssl._create_unverified_context

def gb2faa(phageid, savepath):
    if os.path.exists('%s/%s.faa'%(savepath, phageid)) and os.path.getsize('%s/%s.faa'%(savepath, phageid)) != 0:
        return '%s.faa is existed'%phageid
    savefile = open('%s/%s.faa'%(savepath, phageid), 'w')
    print('now downloading %s.faa'%phageid)
    handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="txt", id=phageid)
    for record in SeqIO.parse(handle, 'gb'):
        count = 0
        for feature in record.features:
            if feature.type == 'CDS':
                count = count + 1
                location = feature.location
                passCharList = ['join', '>', '<']
                checkPassFlag = 0
                for passChar in passCharList:
                    if passChar in str(location):
                        checkPassFlag = 1
                        break
                if checkPassFlag == 1:
                    continue
                if str(location).find('+') != -1:
                    direction = '+'
                elif str(location).find('-') != -1:
                    direction = '-1'
                if feature.type == 'CDS':
                    if 'product' in feature.qualifiers:
                        product = feature.qualifiers['product'][0]
                        if ' ' in product:
                            product = product.replace(' ', '_')
                    else:
                        product = 'unkown'
                    if 'protein_id' in feature.qualifiers:
                        proteinId = feature.qualifiers['protein_id'][0]
                    else:
                        if 'inference' in feature.qualifiers:
                            strInference = str(feature.qualifiers['inference'])
                            if 'RefSeq' in strInference:
                                proteinId = strInference.split('RefSeq:')[1].rstrip(']').rstrip('\'')
                            elif 'SwissProt' in strInference:
                                proteinId = strInference.split('SwissProt:')[1].rstrip(']').rstrip('\'')
                            else:
                                proteinId = 'unknown'
                        else:
                            proteinId = 'unknown'
                    if 'translation' in feature.qualifiers:
                        translation = feature.qualifiers['translation'][0]
                    else:
                        translation = 'unkown'
                    # savefile.write('>' + phageid + '|' + str(proteinId) + '\n')
                    savefile.write('>ref' + '|' + str(proteinId) + '|' + str(count) + '|' + str(location) + '|' + str(product) + '|' + phageid + '\n')
                    if translation[-1] == '\n':
                        savefile.write(translation)
                    else:
                        savefile.write(translation + '\n')

    print('download %s success'%phageid)
    handle.close()