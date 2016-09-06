# This script use otu tables for summaring at given taxonomic level 
# biom convert --to-tsv --header-key taxonomy -i otus.biom -o otus.txt
# Only has been tested using RDP classifier for taxonomy assignment

from optparse import OptionParser
import re, sys
from operator import itemgetter, attrgetter
import numpy as np

def read_input(name):
	metadata = []
	data = []
	file = open(name,"r")
	for line in file:
		if "#" not in line:
			line = line.rstrip().split("\t")
			data.append(line)
		else:
			if "OTU ID" in line.upper():
                                line = re.sub("^#","", line)
				line = line.rstrip().split("\t")
				metadata = line
	file.close()
	return (metadata,data)


def get_taxonomy(data,level):
	RANK = ["k__","p__","c__","o__","f__","g__","s__"]
	pattern = RANK[level-1] + "."
	unclassified = "; ".join(["Unclassified"]*level)
	n_samples = len(data[0])-2
	taxonomy = {}
	
	for line in data:
		tax = line[-1].split("; ")
		if len(tax) >= level:
			if re.search(pattern,tax[level-1]):
				tax ="; ".join(tax[:level])				
				taxonomy[tax] = [list() for i in range(0 , n_samples )]
			else:
				tax = "; ".join([tax[0]] + ["Unclassified"]*(level-1))				
				taxonomy[tax] = [list() for i in range(0 , n_samples )]			
		else:
			tax = "; ".join([tax[0]] + ["Unclassified"]*(level-1))				
			taxonomy[tax] = [list() for i in range(0 , n_samples )]
	
	taxonomy[unclassified] = [list() for i in range(0 , n_samples )]
	
	return taxonomy

def totalReadNumber(taxonomies_otus):
	nsamples = len(taxonomies_otus.itervalues().next())
	total_reads = [0]*nsamples
	
	for key in taxonomies_otus:
		for i in range(0,nsamples):
			total_reads[i] = total_reads[i] + sum(taxonomies_otus[key][i])
	
	return total_reads	

def mymean(otus_nonzero):
	if len(otus_nonzero):
		return round(sum(otus_nonzero) / len(otus_nonzero), 3)
	else:
		return 0

def mypercentage(otus_nonzero, total):
	#return round( (sum(otus_nonzero)/total) * 100 , 3)
	return (sum(otus_nonzero)/total) * 100

def compare_otus(otus):
	nsample = len(otus)
	notus = len(otus[0])
	shared = [1] * notus
	results = [0] * (nsample+1)
		
	for i in range(0,nsample):
		for j in range(0,notus):
			if otus[i][j] > 0:
				results[i]= results[i] + 1
				if shared[j] == 1:
					shared[j] = 1
			else:
				shared[j] = 0
	results[nsample] = shared.count(1)				
	return results	

def analize_otus(taxonomies_otus):
	total_reads = totalReadNumber(taxonomies_otus)
	results = dict()

	for key in taxonomies_otus:
		n_samples = len(taxonomies_otus[key])
		values= [list() for i in range(0,n_samples+1)]
		
		for i in range(0 , n_samples):
			otus = taxonomies_otus[key][i]
			otus_nonzero = filter(lambda a: a != 0, otus)
			
			tax_abundance = sum(otus_nonzero)
			tax_percentage = mypercentage(otus_nonzero, total_reads[i])
			otus_min_count = min(otus_nonzero or [0])
			otus_max_count = max(otus_nonzero or [0])
			otus_mean_count = mymean(otus_nonzero)
			
			values[i] = [ [tax_abundance, tax_percentage], [otus_min_count, otus_max_count, otus_mean_count] ]
		
		values[n_samples] = compare_otus(taxonomies_otus[key])		
		results[key] = values
		
	return results

	
def assig_otus(taxonomies_otus, n_samples, line, tax): # NEW FUNCTION
	for i in range(0,len(line[1:-1])):
		count = float(line[i+1])
		taxonomies_otus[tax][i].append(count)
	return taxonomies_otus

def print_output(taxonomies,metadata): # PARA ELMINAR
	levels=[0,1,5] # show just k,p,g	-> 0:k - 6:s
	sample_names = metadata[1:-1]
	sample_reads = [re.sub("^","#reads ",i) for i in sample_names]
	sample_percentage = [re.sub("^","%",i) for i in sample_names]
	sample_otus = [re.sub("^","#OTUs ",i) for i in sample_names] + ["#OTUs shared"]
	header = ["#Taxonomy"]*len(levels) + sample_reads + sample_percentage + sample_otus
	print "\t".join(header)
	
	data_raw = []
	for key in taxonomies:
		tax_name = np.array(key.split("; "))
		tax_split = list(tax_name[levels])
		data_raw.append( tax_split + taxonomies[key])
		#for join names using "; "
		#tax_split = "; ".join(tax_split)
		#data_raw.append( [tax_split] + taxonomies[key])

	data_sorted = sorted(data_raw, key=itemgetter(0,1,3), reverse=True)
	for line in data_sorted:
		print "\t".join(str(i) for i in line)


def output_format(results, metadata, outputfile, level):
	samples = metadata[1:-1]
	header_tax = ["L"+str(i) for i in range(1,level+1)]

	header_base = ["#sequence_","%sequence_","#OTUs_min_count_","#OTUs_max_count_","#OTUs_mean_count_"]
	header= []
	count=0
	for sample in samples:
		for i in range(0,len(header_base)):
				header = header + [re.sub("$", sample, header_base[i])]
						
	for sample in samples:
		header = header + [re.sub("^","#OTUs_",sample)]
	
	header = header_tax + header + ["#OTUs_shared"]	
	
	file = open(outputfile,"w")
	file.write("\t".join(header) + "\n")
	
	for key in results:
		length = len(results[key])
		string = []
		for i in range(0,length-1):
			sample_length = len(results[key][i])
			for j in range(0, sample_length ):
				string = string + results[key][i][j]
		
		string = key.split("; ") + string + results[key][length-1]
		string = "\t".join([str(x) for x in string])			
		file.write(string + "\n")	
	file.close()

def summarize_taxa(metadata,data,level, outputfile):
	RANK = ["k__","p__","c__","o__","f__","g__","s__"]
	level = level
	pattern = RANK[level-1]+"."
	n_samples = len(metadata[1:-1])
	unclassified = "; ".join(["Unclassified"]*level)
	taxonomies_otus = get_taxonomy(data,level)

	#assign number of reads for taxa with name
	for line in data:
		taxonomy = line[-1].split("; ")

		if len(taxonomy) >= level:
			if re.search(pattern,taxonomy[level-1]):
				tax = "; ".join(taxonomy[:level])
				taxonomies_otus = assig_otus(taxonomies_otus, n_samples, line, tax) #### NEW FUNCTION
				
			else:
				tax = unclassified
				unclass_pattern = RANK[0]+"."			
				if re.search(unclass_pattern, taxonomy[0]):
					tax = taxonomy[0] + "; " +"; ".join(["Unclassified"]*(level-1))
				taxonomies_otus = assig_otus(taxonomies_otus, n_samples, line, tax) #### NEW FUNCTION


		else:
			tax = unclassified
			unclass_pattern = RANK[0]+"."			
			if re.search(unclass_pattern, taxonomy[0]):
				tax = taxonomy[0] + "; " +"; ".join(["Unclassified"]*(level-1))
			taxonomies_otus = assig_otus(taxonomies_otus, n_samples, line, tax) #### NEW FUNCTION


	
	#print_output(taxonomies_otus,metadata)
	results = analize_otus(taxonomies_otus)
	output_format(results, metadata, outputfile, level)


def main():
	parser = OptionParser()
	parser.add_option("-i", "--input-file", dest="inputfile", help="Input file name")
	parser.add_option("-o", "--output-file", dest="outputfile", help="output file name")
	parser.add_option("-L", "--level", dest="level", default=2, type="int", help="select taxonomic level to analize. Default=2 (Phylum)")
	parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True, help="don't print status messages to stdout")

	(options, args) = parser.parse_args()

	if options.inputfile and options.outputfile:
		(metadata, data) = read_input(options.inputfile)
		summarize_taxa(metadata,data,options.level, options.outputfile)
		
	else:
		print "Input file is required"
		sys.exit(0)


if __name__ == "__main__":
	main()
