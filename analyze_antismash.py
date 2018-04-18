
import os, sys
import subprocess
from Bio import SeqIO
import pandas as pd


class SID_records:
	def __init__(self, name):
		self.name = name
		self.genome_file = ""
		self.genbank_file = ""
		self.nrps_recs = []
		self.nrps_CDS = ""
		self.diamond_db = ""
		self.diamond_m8 = ""

def run_antismash(genome_fasta, out_dir_name):
	'''Runs antismash. Takes in genome.fna and name for output folder'''

	subprocess.call(["antismash" , "--outputfolder", out_dir_name, genome_fasta])

def make_fasta_nrps_CDS(SID):
	'''Makes a fasta file when passed a SID object. CDS sequences with translation tags will be pulled out.
	Each header will include the record.id, the feature.type, and the locus_tag.'''
	file = open(SID.name + "_nrps.faa", "w")
	for record in SID.nrps_recs:
		for feature in record.features:
			if "translation" in feature.qualifiers and feature.type == "CDS":
				record_id = record.id
				feature_type = feature.type
				locus_tag = feature.qualifiers["locus_tag"][0]
				translation_seq = feature.qualifiers["translation"][0]
				file.write(">" + record_id + " " + feature_type + " " + locus_tag + "\n")
				file.write(translation_seq + "\n")
	file.close

def run_diamond(SID, query):
	'''Runs diamond using the antismash output files. Pass it an SID object and a query fasta file'''
	for file in os.listdir(start_dir + "/" + SID.name + "/"):
		if ("final.gbk" in file):
			file_location = (start_dir + '/' + SID.name + '/' + file)
			
			SID.genbank = file_location
	
	records = list(SeqIO.parse(SID.genbank, "genbank"))
	for record in records:
		for feature in record.features:
			if feature.type == "cluster" and feature.qualifiers["product"] == ["nrps"]:
				SID.nrps_recs.append(record)
	
	make_fasta_nrps_CDS(SID)
	SID.nrps_CDS = SID.name + "_nrps.faa"
	## Checks to make sure the nrps CDS file isn't empty
	if os.stat(start_dir + "/" + SID.nrps_CDS).st_size != 0:
		subprocess.call(["diamond", "makedb", "--in", SID.nrps_CDS, "-d", SID.name])
		SID.diamond_db = SID.name + ".dmnd"
			
		subprocess.call(["diamond", "blastp", "-d", SID.diamond_db, "-q", query, "-o", start_dir + "/" + SID.name + ".m8", "-f", "6", "qseqid", "sseqid", "qcovhsp", "pident"])
		SID.diamond_m8 = SID.name + ".m8"		

	

# Initializing the start_dir as the directory containing the genome files, file_list as a list
# of all the filenames in that directory, SID_list will contain all the SID objects, the query_fasta
# is used by diamond and needs to be given in command line as the first argument, final_table is a
# dataframe containing the diamond results, QC_cutoff is the query coverage cutoff for hits that will be reported in
# final_table, PI_cutoff is the percent identity cutoff for hits that will be reported in the
# final_table. QC_cutoff and PI_cutoff are given as command line arguments
start_dir = os.getcwd()
file_list = os.listdir(start_dir)
SID_list = []
query_fasta = sys.argv[1]
final_table = pd.DataFrame(columns = ("QSID","SSID","QC","PI"))
QC_cutoff = sys.argv[2]
PI_cutoff = sys.argv[3]
genomes_analyzed = 0

# Iterating though file_list, grabbing all .fna files and making SID objects out of them. Also 
# setting SID.name to the SID of the genome and SID.genome_file to the name of the .fna file. Then
# appending SID_list with the new SID object
for genome in file_list:
	if genome.endswith(".fna"):
		SID_n = genome.split(".")[0]
		SID = SID_records(SID_n)
		SID.genome_file = genome
		SID_list.append(SID)
		


# Running diamond on each SID object
for SID in SID_list:
	run_diamond(SID, query_fasta)
	genomes_analyzed += 1
	print(str(genomes_analyzed) + " genomes analyzed")
	## Checks to make sure the nrps CDS file isn't empty
	if os.stat(start_dir + "/" + SID.nrps_CDS).st_size != 0:
		# Creating a dataframe from the diamond .m8 output file, looking for Percent Identity above certain cutoff
		# and ads them to the final table dataframe
		df = pd.read_csv(start_dir + "/" + SID.diamond_m8, names=["QSID","SSID","QC","PI"], delimiter='\t')
		for i in range(len(df.index - 1)):
			if df.loc[i,"PI"] > int(PI_cutoff) and df.loc[i,"QC"] > int(QC_cutoff):
				final_table = final_table.append(df.loc[i])

# Converting dataframe to .csv file
final_table.to_csv("antismash_results_QC" + QC_cutoff + "_PI" + PI_cutoff + ".csv")
	

	




	
		
