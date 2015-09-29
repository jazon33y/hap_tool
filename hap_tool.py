#===================================================
# This python script was developed by Jason O'Rawe
# as part of his PhD study in 2014 at SBU.
#===================================================

import argparse
from argparse import RawTextHelpFormatter
import os, sys, math, bz2, itertools, re, getopt, numpy
from scipy.stats import beta
import pickle

#======
# Class
#======

class mitchondrial_genome(object):
	def __init__(self, inputfiles = None, type = None, GRCh37 = None, NC_012920p1 = None, hg19 = None, haplogroups = None, poly_hash = None, polys = None):
		self.type = type
		self.files = inputfiles
		self.ref_files = {'GRCh37':GRCh37, 'NC_012920p1':NC_012920p1,'hg19':hg19}
		self.chrM = {}
		self.variants = {}
		self.haplogroups = {}
		self.fasta = {}
		self.expected_variants = {}
		self.observed_variants = {}
		self.false_negative_variants = {}
		self.FN_rate = []
		self.get_chrM()
	def get_chrM(self):
		for i in self.files:
			self.chrM[i]={}
			if i[-3:]=='bz2':
				file = bz2.BZ2File(i,'r')
			else:
				file = open(i,'r')
			if (self.type=='vcf-hg19') or (self.type=='vcf-grch37'):
				col_num = 0
			elif self.type == "var":
				col_num = 3
			ticktock=0
			for line in file.readlines():
				ticktock+=1
				cols = line.split("\t")
				if (line[:1] == "#") or (cols[0]=='\n'):
					continue
				if (cols[col_num][:2] == "MT"):
					cols[col_num] = 'M'
					self.chrM[i][ticktock] = ("\t".join(cols).rstrip('\r\n'))
				if (line.split(str('\t'))[col_num]).rstrip('\r\n')=='M' or (line.split(str('\t'))[col_num]).rstrip('\r\n')=='chrM':
					pass
				else:
					continue
				if (cols[col_num][:3] == "chr"):
					cols[col_num] = cols[col_num][3:] ## removes the characters: 'chr'
					self.chrM[i][ticktock] = ("\t".join(cols).rstrip('\r\n'))
				else:
					self.chrM[i][ticktock] = (line.rstrip('\r\n'))
			file.close()
		return
	def process_chrM(self): 
		chrom_m_old_seq = self.ref_files['GRCh37'] #GRCh37
		chrom_m_seq = self.ref_files['NC_012920p1'] #NC_012920.1
		for outer in self.files:
			expanded_chrM = []
			x = []
			for i in self.chrM[outer].keys():
				x.append(self.chrM[outer][i])
			if self.type=='var':
				for i in x:
					i_s = i.split('\t')
					begin = str(int(i_s[4]) + 1)
					end = str(int(i_s[5]) ) #+1)
					varType = i_s[6]
					reference = chrom_m_seq[int(i_s[4])]
					alleleSeq = i_s[8]

					if i_s[6]=='no-call' or i_s[6]=='no-call-ri' or i_s[6]=='no-call-rc' or i_s[6]=='ref-inconsistent':
						continue
					elif i_s[6]=='snp':
						expanded_chrM.append(begin + '\t' + end + '\t' + varType + '\t' + reference + '\t' + alleleSeq)			
					elif i_s[6]=='sub': # split subs, as they might be SNPs next to each other.  Suggestion by Brenna Henn 2014.
						if len(i_s[7])==len(alleleSeq):
							count = 0
							for ii in range(int(begin),int(end)+1):
								begin = str(ii)
								end = str(ii)
								varType = 'snp'
								reference = chrom_m_seq[ii-1]
								alleleSeq = str(list(i_s[8])[count])
								expanded_chrM.append(begin + '\t' + end + '\t' + varType + '\t' + reference + '\t' + alleleSeq)
								count = count+1			
					elif i_s[6]=='ref':
						start = i_s[4]
						open_end = i_s[5]
						for each in range( int(start),int(open_end) ): #+ 1 ):
							begin = str(each + 1)
							end = str(each + 1)
							varType = 'ref'
							reference = chrom_m_seq[each]
							alleleSeq = chrom_m_old_seq[each]
							expanded_chrM.append(begin + '\t' + end + '\t' + varType + '\t' + reference + '\t' + alleleSeq)  

			if self.type=='vcf-hg19' or self.type=='vcf-grch37':
				if self.type=='vcf-hg19':
					chrom_m_old_seq = self.ref_files['hg19'] #hg19
				list_of_allele_pos = []
				for i in x:
					i_s = i.split('\t')
					if "," in i_s[4]: # split multi allelic sites
						multi_allele = i_s[4].split(",")
						for allele in multi_allele[1:]:
							i_s[4] = allele
							x.append("\t".join(i_s))
						i_s[4] = multi_allele[:1][0]	
					if len(i_s[3]) > 1 or len(i_s[4]) > 1:  # if not a snp
						continue
					elif i_s[3]=='-' or i_s[4]=='-': # if not a snp
						continue
					elif i_s[3]=='.' or i_s[4]=='.': # if not a snp
						continue
					alele_pos = int(i_s[1])
					list_of_allele_pos.append(alele_pos)

					if self.type=='vcf-hg19':
						if alele_pos>=311 and alele_pos<=316:
							alele_pos=alele_pos-1
						elif alele_pos>=318 and alele_pos<=3108:
							alele_pos=alele_pos-2
						elif alele_pos>=3019 and alele_pos<=16190:
							alele_pos=alele_pos-1
						elif alele_pos>=16191 and alele_pos<=16571:
							alele_pos=alele_pos-2
						else:
							alele_pos=alele_pos

					begin = str(int(alele_pos))
					end = begin
					varType = 'snp' #  should be snp at this point
					reference = chrom_m_seq[int(alele_pos)-1] # python is base 0, the above should be base 1
					alleleSeq = i_s[4]
					expanded_chrM.append(begin + '\t' + end + '\t' + varType + '\t' + reference + '\t' + alleleSeq)

				refs = list(set(range(1,len(chrom_m_old_seq)+1)) - set(list_of_allele_pos))
				for i in refs:
					if self.type=='vcf-hg19':
						if i>=311 and i<=316:
							j=i-1
						elif i>=318 and i<=3108:
							j=i-2
						elif i>=3019 and i<=16190:
							j=i-1
						elif i>=16191 and i<=16571:
							j=i-2
						else:
							j=i
					begin = str(i)
					if self.type=='vcf-hg19':
						begin = str(j)
						reference = chrom_m_seq[j-1]
					else:
						begin = str(i)
						reference = chrom_m_seq[i-1] # python is base 0, the above should be base 1
					end = begin
					varType = 'ref' # should be only ref calls now
					alleleSeq = chrom_m_old_seq[i-1] # python is base 0, the above should be base 1
					expanded_chrM.append(begin + '\t' + end + '\t' + varType + '\t' + reference + '\t' + alleleSeq)			
			fixed_chromM = []
			for i in expanded_chrM:
				in_split = i.split('\t')
				if in_split[3]!=in_split[4]:
					fixed_chromM.append(i)

			range_M = []
			variantz = []
			for i in fixed_chromM:
				i_s = i.split('\t')
				i_s_i = int(i_s[0])
				range_M.append(i_s_i)

			for i in fixed_chromM:
				i_s = i.split('\t')
				if ( (i_s[2]=='snp') or (i_s[2]=='ref') ):
					variantz.append(i_s[0] + i_s[4])

			self.variants[outer] = variantz

			def natural_sort(l): #natural sorting algorithm
				convert = lambda text: int(text) if text.isdigit() else text.lower()
				alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
				return sorted(l, key = alphanum_key)

			fasta_data = natural_sort(expanded_chrM)
			fasta = []
			for i in fasta_data:
				fasta.append(i[-1:])
			self.fasta[outer] = '>' + outer + '\n' + ''.join(fasta)

	def get_haplogroups(self):
		for i in self.variants.keys():
			test_list = self.variants[i]
			haplogroup_score = []
			for ii in haplogroups.keys():
				intersect_k = list(set(haplogroups[ii]) & set(test_list))
				weights_k = [poly_hash[iii] for iii in intersect_k]
				sum_weights_k = sum(weights_k) 

				weights_m = [poly_hash[iii] for iii in haplogroups[ii]]
				sum_weights_m = sum(weights_m)

				intersect_n = list(set(polys) & set(test_list))
				weights_n = [poly_hash[iii] for iii in intersect_n]
				sum_weights_n = sum(weights_n) 

				if sum_weights_m==0 or sum_weights_n==0 or sum_weights_k==0:
					continue # sum of the weights cannot be 0, if it is then skip - something went wrong
				score = (float(.5)*(float(sum_weights_k)/sum_weights_m))  + \
						(float(.5)*(float(sum_weights_k)/sum_weights_n) ) # needs to be float value - this could be shorter

				haplogroup_score.append([ii,score])
			max_score = [0,float(0)]
			for ii in range(len(haplogroup_score)):
				test_score = haplogroup_score[ii]
				if test_score[1] > max_score[1]:
					max_score = test_score
			self.haplogroups[i] = max_score

	def get_FN(self):
		wins = []
		losses = []
		all_var_locs = []
		for i in self.variants.keys():
			try:
				a = self.variants[i]
				b = haplogroups[self.haplogroups[i][0]]
				failures = len(set(b) - set(a))
				failed_vars = (set(b) - set(a))
				failed_var_list = list(failed_vars)
				out_list = []
				for ii in range(0,len(failed_var_list)):
					list_of_ele = list(failed_var_list[ii])
					if (("!" in list_of_ele) and (len(list(set(list_of_ele[-2:]))) == len(list_of_ele[-2:]))): # remove back mutation from the expected list
						kill = set([''.join(list_of_ele[:-1])])
						failures_b = len((set(b)- kill) - set(a))
						failed_vars_b = ((set(b)- kill) - set(a))
						failed_var_list_b = list(failed_vars_b)
						b = set(b) - set([''.join(list_of_ele)])
					if ("." in list_of_ele) or ("d" in list_of_ele):
						kill = set([''.join(list_of_ele)])
						failures_b = len((set(b)- kill) - set(a))
						failed_vars_b = ((set(b)- kill) - set(a))
						failed_var_list_b = list(failed_vars_b)
						b = set(b) - set([''.join(list_of_ele)])
				failures = failures_b
				failed_vars = failed_vars_b
				failed_var_list = failed_var_list_b
				for ii in range(0,len(failed_var_list)):
					list_of_ele = list(failed_var_list[ii])
					if (("!" in list_of_ele) and (len(list(set(list_of_ele[-2:]))) == len(list_of_ele[-2:]))): 
						list_of_ele.pop()
						list_of_ele.pop()
						out_list.append("".join(list_of_ele))
					else:
						list_of_ele.pop()
						out_list.append("".join(list_of_ele))
				self.observed_variants[i] = sorted(a)
				self.expected_variants[i] =  sorted(b)
				self.false_negative_variants[i] = failed_var_list
				all_var_locs.extend(out_list)
				sucesses = len(b) - len(out_list)
				failures = len(out_list)
			except:
				print sys.exc_info()
				continue # if there is an error, skip it and go to the next
			wins.append(sucesses)
			losses.append(failures)

		ones = map(int,list(str(1)*sum(losses)))
		zeros = map(int,list(str(0)*sum(wins))) 
		ones.extend(zeros)
		n = 1000
		p1 = 0.5 # Jefferys prior
		p2 = 0.5 # Jefferys prior
		len_ones = len(ones)

		post = []
		for i in range(n):
			post.append(beta.rvs(sum(ones)+p1,len_ones-sum(ones)+p2)) # generate the posterior
		self.FN_rate.append(numpy.mean(post))
		self.FN_rate.extend(beta.ppf([float(0.05),float(.95)],sum(ones)+p1,len_ones-sum(ones)+p2))
	
#======
# Help 
#======

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description="This script uses phylotree information to: \n\n\t"\
													"- return chrM in FASTA format\n\t"\
													"- estimate haplogroup\n\t"\
													"- estimate false negative call rates\n\t"\
													"- return all loci which should have a SNP",formatter_class=RawTextHelpFormatter)

	parser.add_argument('-ana',metavar='<ana>',help='analysis type; fasta, haplogroup, FN_rate, FN_locus', default='haplogroup')
	parser.add_argument('-ft',metavar='<ft>',help='file type; vcf-hg19, vcf-grch37, var', default='vcf-hg19')
	parser.add_argument('-file',metavar='<file>',help='input file, .vcf or .var', nargs='*', required=True)

	args = parser.parse_args()
	ana = args.ana
	ft = args.ft
	file = args.file

	analysis_type = ana
	type = ft
	files = file

	#========================
	# Load *required* objects
	# 
	# Note that these are not
	# necessarily computed 
	# beforehand, this just makes
	# the computation time very fast
	# 
	#
	obj1 = open("objs/haplogroups.obj",'rb')
	obj2 = open("objs/poly_hash.obj",'rb')
	obj3 = open("objs/polys.obj",'rb')
	obj4 = open("objs/GRCh37.obj",'rb')
	obj5 = open("objs/NC_012920p1.obj",'rb')
	obj6 = open("objs/hg19.obj",'rb')
	# 
	# 
	haplogroups = pickle.load(obj1)
	poly_hash = pickle.load(obj2)
	polys = pickle.load(obj3)
	GRCh37 = pickle.load(obj4)
	NC_012920p1 = pickle.load(obj5)
	hg19 = pickle.load(obj6)
	# 
	# 
	obj1.close()
	obj2.close()
	obj3.close()
	obj4.close()
	obj5.close()
	obj6.close()
	# 
	#========================

	#========
	# do work
	#========

	mt_Genomes = mitchondrial_genome(inputfiles = files, type = type, GRCh37 = GRCh37 , NC_012920p1 = NC_012920p1, hg19 = hg19,haplogroups = haplogroups, poly_hash = poly_hash, polys = polys)

	mt_Genomes.process_chrM()
	mt_Genomes.get_haplogroups()
	mt_Genomes.get_FN()

	if analysis_type == "fasta":
		for i in mt_Genomes.fasta.keys():
			print mt_Genomes.fasta[i]

	elif analysis_type == "haplogroup":
		for i in mt_Genomes.haplogroups.keys():
			print str(i) + ": " + str(mt_Genomes.haplogroups[i][0]) + "; " + str(mt_Genomes.haplogroups[i][1])

	elif analysis_type == "FN_rate":
		print "\n\'Error rate:\'"
		print "Expectation: ",mt_Genomes.FN_rate[0]
		print "95% credible interval: [",mt_Genomes.FN_rate[1],",",mt_Genomes.FN_rate[2],"]\n"

	elif analysis_type == "FN_locus":
		for i in mt_Genomes.observed_variants.keys():
			print 'Locations where calls should be made but were not, ' + i + ':'#,map(int,all_var_locs),'\n'
			print 'observed ',i,':\t', mt_Genomes.observed_variants[i]
			print 'expected ',i,':\t', mt_Genomes.expected_variants[i]
			print 'False negatives:\t',mt_Genomes.false_negative_variants[i]

	else: print "\nInvalid analysis type\n"

