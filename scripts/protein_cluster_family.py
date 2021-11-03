import os
import xml.etree.ElementTree as ET
import numpy as np
import threading
import argparse
import re

def mkdir(dirname):
	command = "mkdir -p "+dirname
	os.system(command)

def MMseqs2_cluster_easy(fasta_file,outfile,temp_dir,sensitivity,identity=0.5,coverage=0.8,evalue=1e-3):
	command  = "mmseqs easy-cluster %s %s %s -s %s --min-seq-id %s -c %s -e %s --cluster-mode 0"%(fasta_file,outfile,temp_dir,str(sensitivity),str(identity),str(coverage),str(evalue))
	os.system(command)
	print(command)

def MMseqs2_search_two_easy(query_fasta_file,target_fasta_file,outfile,temp_dir,sensitivity=0.7,identity=0.5,coverage=0.5,evalue=1e-3):
	command = "mmseqs easy-search %s %s %s %s -s %s --min-seq-id %s -c %s --e-profile %s"%(query_fasta_file,
		target_fasta_file,outfile,temp_dir,
		str(sensitivity),str(identity),str(coverage),str(evalue))
	os.system(command)
	print(command)

def mmseqs_create_db(query_fasta_file,target_DB,temp_dir):
	command = "mmseqs createdb %s %s"%(query_fasta_file,target_DB)
	os.system(command)
	command = "mmseqs createindex %s %s"%(target_DB,temp_dir)
	os.system(command)

def mmseqs_cluster(DB,outfile,temp_dir,sensitivity=7.5,identity=0.4,coverage=0.7,evalue=1e-3):
	command = "mmseqs cluster %s %s %s -s %s --min-seq-id %s -c %s -e %s --cluster-mode 0"%(DB,outfile,temp_dir,str(sensitivity),str(identity),str(coverage),str(evalue))
	os.system(command)

def mmseqs_msa(DB,cluster_file,msa_file):	
	command = "mmseqs result2msa %s %s %s %s"%(DB,DB,cluster_file,msa_file)
	os.system(command)

def build_hmm_database_from_msa(msa_file,DB,HHsuite_db):
	compile_dir = "/100T/software/hhsuite/hh-suite/build/compile/data"
	command1 = "ln -s %s_h %s_header.ffdata"%(DB,HHsuite_db)
	command2 = "ln -s %s_h.index %s_header.ffindex"%(DB,HHsuite_db)
	command3 = "ln -s %s %s_sequence.ffdata"%(DB,HHsuite_db)
	command4 = "ln -s %s.index %s_sequence.ffindex"%(DB,HHsuite_db)
	#command5 = "mpirun -np 2 cstranslate_mpi -i %s -o %s/%s -A %s/cs219.lib -D %s/context_data.lib -x 0.3 -c 4 -I ca3m -b"%(msa_file,HHsuite_db_dir,prefix,compile_dir,compile_dir)
	command5 = "cstranslate -i %s -o %s -A %s/cs219.lib -D %s/context_data.lib -x 0.3 -c 4 -I ca3m -b"%(msa_file,HHsuite_db,compile_dir,compile_dir)
	
	os.system(command1)
	os.system(command2)
	os.system(command3)
	os.system(command4)
	os.system(command5)
	print(command5)

def hhblits_compare(msa_file,outfile,a3m_file,HHsuite_db,cpu_num=10):
	command = "hhblits -i %s -o %s -oa3m %s -d %s -v 0 -p 50 -z 4 -Z 32000 -B 0 -b 0 -cpu %s"%(msa_file,outfile,a3m_file,HHsuite_db,cpu_num)
	os.system(command)

def mafft(query_file,outfile,thread_num):
	command = "mafft --thread %s %s > %s"%(thread_num,query_file,outfile)
	os.system(command)

def build_hmm_hmmbuild(msa_file,hmmfile,hmm_name):
	command = "hmmbuild -n %s %s %s"%(hmm_name,hmmfile,msa_file)
	os.system(command)

def build_hmm_hmmpress(hmm_file):	
	command = "hmmpress %s"%hmm_file
	os.system(command)

class MyThread_hmmbuild(threading.Thread):
	def __init__(self,msa_file,hmmfile,hmm_name):
		threading.Thread.__init__(self)
		self.msa_file = msa_file
		self.hmmfile = hmmfile
		self.hmm_name = hmm_name
	def run(self):
		build_hmm_hmmbuild(self.msa_file,self.hmmfile,self.hmm_name)

def rm_h3_file(hmm_file):
	file_list = ['h3p','h3m','h3i','h3f']
	root_path = '/'.join(hmm_file.split('/')[0:-1])
	for file_item in file_list:
		cur_file = '%s.%s'%(hmm_file,file_item)
		if os.path.exists(cur_file):
			cmd_rm_h3_file = 'rm -r %s/*.h3*'%(root_path)
			os.system(cmd_rm_h3_file)

def build_cluster_hmm(msa_file,outdir,save_prefix,thread_num):
	save_summary_file = save_prefix+'_summary.txt'
	save_represent_file = save_prefix+'_ref.faa'
	save_hmm_file = save_prefix+'_hmm'
	f_save = open(save_summary_file,'w')
	f_save.write('cluster_id\tcluster_size\trepresent_id\tgroup_id\n')
	f_save_ref = open(save_represent_file,'w')
	f_save_hmm = open(save_hmm_file,'w')
	with open(msa_file) as f:
		contents = f.read().split('\n\x00>')
	#print(contents)
	cluster_id = 0
	#print(len(contents))
	tsk = []	
	for cluster in contents:
		cluster_id = cluster_id+1
		cluster_size = len(cluster.strip().strip('>').split('>'))
		c_dir = os.path.join(outdir,'cluster'+str(cluster_id))
		mkdir(c_dir)
		c_msa_file = os.path.join(c_dir,'cluster'+str(cluster_id)+'_msa')
		with open(c_msa_file,'w') as f:
			f.write('>'+cluster.strip().strip('>').strip()+'\n')
		c_hmm_file = os.path.join(c_dir,'cluster'+str(cluster_id)+'_hmm')
		if cluster_size>1:
			try:
				tsk.append(MyThread_hmmbuild(c_msa_file,c_hmm_file,cluster_id))
			except:
				print('Error: unable to start thread')
	for t in tsk:
		t.start()
		while True:
			if (len(threading.enumerate()) <= int(thread_num)):
				break
	for t in tsk:
		t.join()
	
	cluster_id = 0
	#print(len(contents))	
	for cluster in contents:
		cluster_id = cluster_id+1
		cluster_size = len(cluster.strip('>').split('>'))
		c_dir = os.path.join(outdir,'cluster'+str(cluster_id))
		mkdir(c_dir)
		c_msa_file = os.path.join(c_dir,'cluster'+str(cluster_id)+'_msa')
		with open(c_msa_file,'w') as f:
			f.write('>'+cluster.strip().strip('>').strip()+'\n')
		c_hmm_file = os.path.join(c_dir,'cluster'+str(cluster_id)+'_hmm')
		if cluster_size>1:
			with open(c_hmm_file) as f:
				f_save_hmm.write(f.read().strip()+'\n')
				f_save_hmm.flush()
		index = 0
		for protein in cluster.strip('>').split('>'):
			protein_id = protein.strip().split('\n')[0].strip().strip('>').strip()
			if index==0:
				represent_id = protein_id
				index = index+1
			f_save.write(str(cluster_id)+'\t'+str(cluster_size)+'\t'+represent_id+'\t'+protein_id+'\n')
			f_save.flush()
			if cluster_size>1:
				f_save_ref.write('>'+str(cluster_id)+'|'+protein.strip().replace('-','')+'\n')
				f_save_ref.flush()
	
	f_save.close()
	f_save_ref.close()
	f_save_hmm.close()
	rm_h3_file(save_hmm_file)
	build_hmm_hmmpress(save_hmm_file)

def rpsblastp11proteins(fileName,outFileName,domain_file):
	blast_cline = "/100T/pipeline/crispr/tools/blast-2.2.31-1/bin/rpsblast -query %s -comp_based_stats 0 -max_target_seqs 1 -evalue 0.01 -seg no -outfmt 5 -num_threads 20 -db /100T/pipeline/crispr/db/cdd/Cdd -out %s" % (fileName,outFileName)
	os.system(blast_cline)

	text = open(outFileName).read()
	text = re.sub(u"[\x00-\x08\x0b-\x0c\x0e-\x1f]+", u" ", text)
	root = ET.fromstring(text)
	BlastOutput_iterations = root.find("BlastOutput_iterations")
	f = open(domain_file,'w')
	for Iteration in BlastOutput_iterations.findall("Iteration"):
		strTemp = str(Iteration.find("Iteration_query-def").text)
		Iteration_hits = Iteration.find("Iteration_hits")
		for Hit in Iteration_hits.findall("Hit"):
			strDef = str(Hit.find("Hit_def").text)
			Hit_len =  Hit.find("Hit_len").text
			Hit_ID = Hit.find("Hit_id").text
			Hsp = Hit.find("Hit_hsps").find("Hsp")
			Hit_from = Hsp.find("Hsp_hit-from").text
			Hit_to = Hsp.find("Hsp_hit-to").text
			Hsp_evalue = Hsp.find("Hsp_evalue").text
			Hsp_score = Hsp.find("Hsp_score").text
			hit_length = int(Hit_to)-int(Hit_from)+1
			coverage = str(float(hit_length)/ float(Hit_len))
			identity = str(Hsp.find("Hsp_identity").text)
			f.write(strTemp+'\t'+Hit_ID+'\t'+strDef+'\t'+Hsp_evalue+'\t'+identity+'\t'+coverage+'\n')
			f.flush()
	f.close()

def hmmscan_fasta_to_hmm(query_fasta_file,hmmfile,hit_file,hit_parse_file,cluster_id):
	cmd = 'hmmscan -E 1e-5 --tblout %s %s %s > /dev/null' % (hit_file,hmm_file,query_fasta_file)
	os.system(cmd)
	f_save = open(hit_parse_file,'w')
	with open(hit_file) as f:
		contents = f.readlines()
	for line in contents:
		if line[0] == '#':
			continue
		linetab = line.split()
		query_pro = linetab[2]
		#query_cluster_id = query_pro.split('|')[0]		
		hit_cluster_id = linetab[0]
		src_evalue = linetab[4]
		src_score = linetab[5]
		f_save.write(cluster_id+'\t'+hit_cluster_id+'\t'+query_pro+'\t'+src_evalue+'\t'+src_score+'\n')
		f_save.flush()
	f_save.close()

class MyThread_hmmscan(threading.Thread):
	def __init__(self,query_fasta_file,hmmfile,hit_file,hit_parse_file,cluster_id):
		threading.Thread.__init__(self)
		self.query_fasta_file = query_fasta_file
		self.hmmfile = hmmfile
		self.hit_file = hit_file
		self.hit_parse_file = hit_parse_file
		self.cluster_id = cluster_id	
	def run(self):
		hmmscan_fasta_to_hmm(self.query_fasta_file,
			self.hmmfile,self.hit_file,
			self.hit_parse_file,self.cluster_id)

def hmmsearch_hmm_to_fasta(query_hmm_file,fasta_file,hit_file,hit_parse_file):
	cmd = "hmmsearch -E 1e-5 --tblout %s %s %s > /dev/null"%(hit_file,query_hmm_file, fasta_file)
	os.system(cmd)
	f_save = open(hit_parse_file,'w')
	with open(hit_file) as f:
		contents = f.readlines()
	for line in contents:
		if line[0] == '#':
			continue
		linetab = line.split()
		query_hmm = linetab[2]
		hit_pro = linetab[0]
		hit_cluster_id = hit_pro.split('|')[0]		
		src_evalue = float(genequery[linetab[2]][2])
		src_score = float(genequery[linetab[2]][3])
		src_bias = float(genequery[linetab[2]][4])
		dst_evalue = float(linetab[4])
		dst_score = float(linetab[5])
		dst_bias = float(linetab[6])		
		f_save.write(query_hmm+'\t'+hit_cluster_id+'\t'+hit_pro+'\t'+src_evalue+'\t'+src_score+'\t'+dst_evalue+'\t'+dst_score+'\n')
		f_save.flush()
	f_save.close()

def get_mcl_abc_file(hmm_dir,hmm_file,cluster_info_file,matrix_abc_file,hmmscan_dir,thread_num):
	f_save = open(matrix_abc_file,'w')
	with open(cluster_info_file) as f:
		contents = f.readlines()

	tsk = []
	cluster_ids = []
	for line in contents[1:]:
		line = line.strip().split('\t')
		cluster_id = line[0]
		cluster_size = line[1]
		if int(cluster_size)>1:
			#c_hmm_file = os.path.join(msa_dir,'cluster'+cluster_id,'cluster'+cluster_id+'_hmm')
			c_msa_file = os.path.join(hmm_dir,'cluster'+cluster_id,'cluster'+cluster_id+'_msa')
			c_hmmscan_dir = os.path.join(hmmscan_dir,'cluster'+cluster_id)
			mkdir(c_hmmscan_dir)
			c_hit_file = os.path.join(c_hmmscan_dir,'cluster'+cluster_id+'_hmmsearch')
			c_hit_parse_file = os.path.join(c_hmmscan_dir,'cluster'+cluster_id+'_hmmsearch_hit.txt')				
			# hmmsearch_hmm_to_fasta(c_hmm_file,ref_pro_file,c_hit_file,c_hit_parse_file)
			# hmmscan_fasta_to_hmm(c_msa_file,hmm_file,c_hit_file,c_hit_parse_file)
			if cluster_id not in cluster_ids:
				cluster_ids.append(cluster_id)
				try:
					tsk.append(MyThread_hmmscan(c_msa_file,hmm_file,c_hit_file,c_hit_parse_file,cluster_id))
				except:
					print('Error: unable to start thread')
	
	for t in tsk:
		t.start()
		while True:
			if (len(threading.enumerate()) <= int(thread_num)):
				break
	for t in tsk:
		t.join()

	for line in contents[1:]:
		line = line.strip().split('\t')
		cluster_id = line[0]
		cluster_size = line[1]
		if int(cluster_size)>1:
			#c_hmm_file = os.path.join(msa_dir,'cluster'+cluster_id,'cluster'+cluster_id+'_hmm')
			c_msa_file = os.path.join(hmm_dir,'cluster'+cluster_id,'cluster'+cluster_id+'_msa')
			c_hmmscan_dir = os.path.join(hmmscan_dir,'cluster'+cluster_id)
			mkdir(c_hmmscan_dir)
			c_hit_file = os.path.join(c_hmmscan_dir,'cluster'+cluster_id+'_hmmsearch')
			c_hit_parse_file = os.path.join(c_hmmscan_dir,'cluster'+cluster_id+'_hmmsearch_hit.txt')				

			c_cluster_info_dict = {}
			with open(c_hit_parse_file) as f:
				hit_contents = f.readlines()
			for line1 in hit_contents[1:]:
				line1 = line1.strip().split('\t')
				query_pro_id = line1[2]
				hit_cluster_id = line1[1]
				score = line1[4]
				evalue = line1[3]
				if hit_cluster_id not in c_cluster_info_dict.keys():
					c_cluster_info_dict.update({hit_cluster_id:[]})
				c_cluster_info_dict[hit_cluster_id].append([query_pro_id,score,evalue])
			for hit_cluster_id in c_cluster_info_dict.keys():
				hit_pros = c_cluster_info_dict[hit_cluster_id]
				hit_pro_num = len(hit_pros)
				mean_score = np.mean([float(item[1]) for item in hit_pros])
				f_save.write(cluster_id+'\t'+hit_cluster_id+'\t'+str(mean_score)+'\n')
				f_save.flush()
	f_save.close()

def mcl(matrix_abc_file,outdir,prefix):
	tab_file = os.path.join(outdir,prefix+'_tab.txt')
	mci_file = os.path.join(outdir,prefix+'_mci.txt')
	command1 = "mcxload -abc %s --stream-mirror -write-tab %s -o %s"%(matrix_abc_file,tab_file,mci_file)
	os.system(command1)
	print(command1)
	
	mcl_outfile = os.path.join(outdir,prefix+'_mcl.txt')
	command2 = "mcl %s --abc -I 2 -o %s"%(matrix_abc_file,mcl_outfile)
	os.system(command2)
	print(command2)
	
	# outfile = os.path.join(outdir,prefix+'_mcl_result.txt')
	# command3 = "mcxdump -icl %s -tabr %s -o %s"%(mcl_outfile,tab_file,outfile)
	# os.system(command3)
	# print(command3)

def get_protein_family(ori_table_file,mcl_file,cluster_info_file,save_file):
	f_save = open(save_file,'w')
	f_save.write('family_id\tfamily_member\tcluster_id\tcluster_flag\tprotein_id\tprotein_size\tassembly_id\tgenome_id\tgenome_def\tdbscan-swa_range\tphagefinder_range\tmerge_range\tcrispr_array_location_merge\tcrispr_pred_method\tarray_in_prot\tprot_within_array_20000\tcrispr_type_by_cas_prot\tconsensus_repeat\tspacer_max_num\tarray_status\t\n')
	
	table_info_dict = {}
	with open(ori_table_file) as f:
		contents = f.readlines()
	header = contents[0].strip('\n').split('\t')
	for line in contents[1:]:
		line = line.strip().split('\t')
		# print(len(line))
		protein_id = line[0].split('|')[0]
		protein_size = line[1]
		phage_id = line[header.index('phage_id')]
		prot_def = line[header.index('protein_def')]
		prot_start_end = line[header.index('protein_start_end')]
		phage_family = line[header.index('phage_family')]
		prot_index = line[header.index('count')]
		# dbscan_range = line[header.index('dbscan-swa_range')]
		# phagefinder_ranger = line[header.index('phagefinder_range')]
		# merge_range = line[header.index('merge_range')]
		# crispr_array_location_merge = line[header.index('crispr_array_location_merge')]
		# crispr_pred_method = line[header.index('crispr_pred_method')]
		# array_in_prot = line[header.index('array_in_prot')]
		# prot_within_array_20000 = line[header.index('prot_within_array_20000')]
		# crispr_type_by_cas_prot = line[header.index('crispr_type_by_cas_prot')]
		# consensus_repeat = line[header.index('consensus_repeat')]
		# spacer_locus_num = line[header.index('spacer_locus_num')]
		# spacer_max_num = line[header.index('spacer_max_num')]
		# array_status = line[header.index('array_status')]
		if protein_id not in table_info_dict.keys():
			table_info_dict.update({protein_id: [protein_size,
												 '', '', '', '', '']})
		table_info_dict[protein_id][1] = table_info_dict[protein_id][1] + phage_id
		table_info_dict[protein_id][2] = table_info_dict[protein_id][2] + prot_def
		table_info_dict[protein_id][3] = table_info_dict[protein_id][3] + prot_start_end
		table_info_dict[protein_id][4] = table_info_dict[protein_id][4] + phage_family
		table_info_dict[protein_id][5] = table_info_dict[protein_id][5] + prot_index
		# table_info_dict[protein_id][6] = table_info_dict[protein_id][6] + merge_range + '|'
		# table_info_dict[protein_id][7] = table_info_dict[protein_id][7] + crispr_array_location_merge + '|'
		# table_info_dict[protein_id][8] = table_info_dict[protein_id][8] + crispr_pred_method + '|'
		# table_info_dict[protein_id][9] = table_info_dict[protein_id][9] + array_in_prot + '|'
		# table_info_dict[protein_id][10] = table_info_dict[protein_id][10] + prot_within_array_20000 + '|'
		# table_info_dict[protein_id][11] = table_info_dict[protein_id][11] + crispr_type_by_cas_prot + '|'
		# table_info_dict[protein_id][12] = table_info_dict[protein_id][12] + consensus_repeat + '|'
		# table_info_dict[protein_id][13] = table_info_dict[protein_id][13] + spacer_locus_num + '|'
		# table_info_dict[protein_id][14] = table_info_dict[protein_id][14] + spacer_max_num + '|'
		# table_info_dict[protein_id][15] = table_info_dict[protein_id][15] + array_status + '|'
	cluster_info_dict = {}
	with open(cluster_info_file) as f:
		contents = f.readlines()
	for line in contents[1:]:
		line = line.strip().split('\t')
		cluster_id = line[0]
		represent_id = line[2]
		group_id = line[3]	
		if cluster_id not in cluster_info_dict.keys():
			cluster_info_dict.update({cluster_id:[]})
		if represent_id not in cluster_info_dict[cluster_id]:
			cluster_info_dict[cluster_id].append(represent_id)
		if group_id not in cluster_info_dict[cluster_id]:
			cluster_info_dict[cluster_id].append(group_id)
	
	mcl_info_dict = {}
	in_family_proteins = []
	with open(mcl_file) as f:
		contents = f.readlines()
	family_id = 0
	for line in contents:
		c_family = line.strip().split('\t')
		family_id = family_id+1
		c_family_number = [len(list(cluster_info_dict[cluster_id])) for cluster_id in c_family]
		c_family_number = np.sum(c_family_number)
		for cluster_id in c_family:
			if len(list(cluster_info_dict[cluster_id]))>1:
				cluster_flag = 'cluster'
			else:
				cluster_flag = 'singleton'
			
			for group_protein in cluster_info_dict[cluster_id]:
				protein_id = group_protein.split('|')[0]
				in_family_proteins.append(group_protein)
				#mcl_info_dict.update({protein_id:str(family_id)+'\t'+str(c_family_number)+'\t'+cluster_id+'\t'+cluster_flag+'\t'+group_protein+'\t'+'\t'.join(table_info_dict[protein_id]).strip()})
				f_save.write(str(family_id)+'\t'+str(c_family_number)+'\t'+cluster_id+'\t'+cluster_flag+'\t'+group_protein+'\t'+'\t'.join(table_info_dict[protein_id]).strip()+'\n')
				f_save.flush()
				in_family_proteins.append(protein_id)
	
	for cluster_id in cluster_info_dict.keys():
		family_id = family_id+1
		if len(list(cluster_info_dict[cluster_id]))>1:
			cluster_flag = 'cluster'
		else:
			cluster_flag = 'singleton'
		for group_protein in cluster_info_dict[cluster_id]:
			protein_id = group_protein.split('|')[0]
			if group_protein not in in_family_proteins:				
				f_save.write(str(family_id)+'\t'+str(1)+'\t'+cluster_id+'\t'+cluster_flag+'\t'+group_protein+'\t'+'\t'.join(table_info_dict[protein_id]).strip()+'\n')
				f_save.flush()
	f_save.close()

def rm_dir(dir_path):
	cmd_rm = 'rm -rf %s'%dir_path
	os.system(cmd_rm)

if __name__=='__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', help='absPath of genome file.\n')
	parser.add_argument('--output', help='the out dir.\n')
	parser.add_argument('--table', help='protein_info file')
	parser.add_argument('--prefix', help='save prefix')
	args = parser.parse_args()
	if args.input:
		query_fasta_file = args.input
	if args.output:
		outdir = args.output
	if args.table:
		table_file = args.table
	if args.prefix:
		prefix = args.prefix
	else:
		prefix = os.path.basename(outdir)
	
	#step1:mmseqs to cluster proteins and create MSA files
	mmseqs_dir = os.path.join(outdir,'mmseqs')	
	mkdir(mmseqs_dir)
	temp_dir = os.path.join(mmseqs_dir,'temp')
	mkdir(temp_dir)
	
	##create mmseqs db for query fasta
	DB_dir = os.path.join(mmseqs_dir,'DB')
	mkdir(DB_dir)
	DB_file = os.path.join(DB_dir,prefix+'_mmseqs_db')
	mmseqs_create_db(query_fasta_file,DB_file,temp_dir)
	
	##mmseqs cluster
	cluster_dir = os.path.join(mmseqs_dir,'cluster')
	if os.path.exists(cluster_dir):
		rm_dir(cluster_dir)
	mkdir(cluster_dir)
	cluster_file = os.path.join(cluster_dir,prefix+'_cluster')
	mmseqs_cluster(DB_file,cluster_file,temp_dir,7.5,0,0.5,1e-3)
	
	##multiple sequence alignment
	msa_dir = os.path.join(mmseqs_dir,'msa')
	mkdir(msa_dir)
	msa_file = os.path.join(msa_dir,prefix+'_msa')
	mmseqs_msa(DB_file,cluster_file,msa_file)

	#step2:build HMM profiles based on MSA results
	save_cluster_summary_file = os.path.join(mmseqs_dir,prefix+'_summary.txt')
	save_hmm_prefix = os.path.join(mmseqs_dir,prefix)
	hmm_dir = os.path.join(outdir,'hmm')
	mkdir(hmm_dir)
	build_cluster_hmm(msa_file,hmm_dir,save_hmm_prefix,5)

	#step3:hmmscan between each cluster to HMM profiles and get the weight score for MCL
	hmm_file = save_hmm_prefix+'_hmm'
	cluster_info_file = save_hmm_prefix+'_summary.txt'
	mcl_dir = os.path.join(outdir,'mcl')
	mkdir(mcl_dir)
	matrix_abc_file = os.path.join(mcl_dir,prefix+'_mcl_abc.txt')
	hmmscan_dir = os.path.join(outdir,'hmmscan')
	mkdir(hmmscan_dir)
	get_mcl_abc_file(hmm_dir,hmm_file,save_cluster_summary_file,matrix_abc_file,hmmscan_dir,5)

	#step4:MCL to clsuter subclusters into protein families 
	mcl(matrix_abc_file,mcl_dir,prefix)

	#step5:parse results and get result tables
	mcl_file = os.path.join(mcl_dir,prefix+'_mcl.txt')
	# save_file = os.path.join(mcl_dir,prefix+'_result.txt')
	save_file = '%s/%s_result.txt'%(outdir,prefix)
	get_protein_family(table_file,mcl_file,save_cluster_summary_file,save_file)



