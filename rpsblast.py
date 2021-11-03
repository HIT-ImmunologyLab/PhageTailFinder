import os
import xml.etree.ElementTree as ET
import re
import threading

def ParseResult(inFileName,outFileName):
    fout = open(outFileName,'w')
    headers = "uniprot_prot_id\tuniprot_prot_name\tuniprot_prot_def\tpfam_id\tdomain_id\tdomain_name\tdomain_def\tquery_from\tquery_to\thit_from\thit_to\tidentity\te-value\tdomain_coverage\n"
    fout.write(headers)
    fout.flush()

    # 过滤xml文件中的非法字符
    text = open(inFileName).read()
    text = re.sub(u"[\x00-\x08\x0b-\x0c\x0e-\x1f]+", u" ", text)
    root = ET.fromstring(text)

    BlastOutput_iterations = root.find("BlastOutput_iterations")

    for Iteration in BlastOutput_iterations.findall("Iteration"):
        query_info = str(Iteration.find("Iteration_query-def").text)
        uniprot_prot_id = query_info.split(' ')[0]
        uniprot_prot_name = ' '.join(query_info.split('[')[0].split(' ')[1:])
        uniprot_prot_phage =  query_info.split('[')[1].split(']')[0]
        Iteration_query_len = Iteration.find('Iteration_query-len').text
        Iteration_hits = Iteration.find("Iteration_hits")
        for Hit in Iteration_hits.findall("Hit"):
            Hit_length = Hit.find("Hit_len").text
            Hsp = Hit.find("Hit_hsps").find("Hsp")
            Hsp_evalue = Hsp.find("Hsp_evalue").text
            query_from = Hsp.find("Hsp_query-from").text
            query_to = Hsp.find("Hsp_query-to").text
            Hit_from = Hsp.find("Hsp_hit-from").text
            Hit_to = Hsp.find("Hsp_hit-to").text
            identity = str(float(Hsp.find("Hsp_identity").text) / float(Hsp.find("Hsp_align-len").text))
            # coverage = str(float(Hsp.find('Hsp_align-len').text) / float(Iteration_query_len))
            Hit_align_length = abs(float(Hit_to)-float(Hit_from))+1
            coverage = str(float(Hit_align_length) / float(Hit_length))
            strHitID = str(Hit.find("Hit_id").text)
            Domain_ID = strHitID.split('|')[2]
            strHitDef = str(Hit.find("Hit_def").text)
            Domain_Name = strHitDef.split(',')[0]
            Domain_Def = strHitDef.strip(Domain_Name).strip(', ').split('.')[0]

            strWrite = uniprot_prot_id+'\t'+uniprot_prot_name+'\t'+uniprot_prot_phage+'\t'+Domain_ID+'\t'+Domain_Name+'\t'+Domain_Def+'\t'+query_from+'\t'+query_to+'\t'+Hit_from+'\t'+Hit_to+'\t'+identity+'\t'+Hsp_evalue+'\t'+coverage+'\n'
            fout.write(strWrite)
            fout.flush()
    fout.close()

def rpsblast(queryFileName, rpsblastOut):
    xml_out_file = '%s.xml'%rpsblastOut
    blast_cline = "rpsblast -query %s -comp_based_stats 0 -evalue 0.01 -seg no -outfmt 5 -num_threads 20 -db /home/data/CDD_2020_07_22/Cdd -out %s" % (queryFileName, xml_out_file)
    os.system(blast_cline)
    ParseResult(xml_out_file, rpsblastOut+'_rpsblastresult')

class MyRpsThread(threading.Thread):
    def __init__(self,queryFileName, rpsblastOut):
        threading.Thread.__init__(self)
        self.queryFileName = queryFileName
        self.rpsblastOut = rpsblastOut
    def run(self):
        rpsblast(self.queryFileName,self.rpsblastOut)

def runrpsblast(faa_path, rpsblastout_path):

    print('-----------------------------------------------------------------------\n')
    print('--------------------------NOW RUNNING rpsblast-------------------------\n')
    print('-----------------------------------------------------------------------\n')
    os.system('mkdir -p %s' % rpsblastout_path)
    tsk = []
    for root, subdirs, files1 in os.walk(faa_path):
        for subdir in subdirs:
            subdir_path = os.path.join(root, subdir)
            for root1, subdirs1, files in os.walk(subdir_path):
                for file in files:
                    faa_file = os.path.join(root1, file)
                    faa_id = file.split('.')[0]
                    rpsout = rpsblastout_path + '/' + faa_id
                    if os.path.exists(rpsout):
                        pass
                    else:
                        try:
                            tsk.append(MyRpsThread(faa_file, rpsout))
                        except:
                            print('Error: unable to start thread')
    for t in tsk:
        t.start()
        while True:
            if (len(threading.enumerate()) <= int(8)):
                break
    for t in tsk:
        t.join()
