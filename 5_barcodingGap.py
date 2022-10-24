#!/usr/bin/env python

"""

NAME
	barcoding_gap.py

SYNOPSIS
	barcoding_gap.py -f fasFile -t taxFile

PARAMETERS
	-f
		Fasta file with accession numbers (.fasta)
	-t
		Taxonomy file with the equivalences between accession numbers and taxonomy (superphylum;phylum;class;order;family;genus;species)
"""
def  main():

#		import numpy
	import sys
	import glob, os
#        import re
#        import gzip
#        from Bio.Data.IUPACData import ambiguous_dna_values
#        from Bio.Seq import Seq

# Test if the number of arguments is rigth

	if len(sys.argv) < 2:
		print __doc__
		sys.exit(1)

# Process arguments

	fasFile = ""
	taxFile = ""
	argv = sys.argv[1:]

	while len(argv) > 0 :
        	arg = argv[0]
        	if arg == "-f" :
        		if len(argv) > 1 :
            			fasFile = argv[1]
	            	if os.path.exists(fasFile) == 0 :
        	    		print("\n\tERROR: Fasta file \"" + fasFile + "\" could not be found!\n")
                		sys.exit(0)
            		argv = argv[2:]
	            	continue

		if arg == "-t" :
    			if len(argv) > 1 :
            			taxFile = argv[1]
	            	if os.path.exists(taxFile) == 0 :
        	    		print("\n\tERROR: Taxonomy file \"" + taxFile + "\" could not be found!\n")
                		sys.exit(0)
            		argv = argv[2:]
	            	continue
# calculate distances among sequences

	disFile = fasFile.replace(".fasta", ".dist")	
	print "Blast output will be printed to " + disFile + "\n"
	cmd = "makeblastdb -in  " + fasFile + " -dbtype nucl"
	print cmd
    	os.system(cmd)
	cmd = "blastn -db " + fasFile + " -outfmt '6 qseqid sseqid pident' -out " + disFile + " -query " + fasFile + " -num_threads 8"
	print cmd
        os.system(cmd)
	print "Blast process finished\n"

# create taxonomy dictionary

	taxLines = open(taxFile, "r").read().splitlines()
	accTax = {}
	for line in taxLines:
		accTax[line.split("\t")[0]] = line.split("\t")[1]


# count number of comparisons

        
	print "Counting pairs" + taxFile + "\n"
	i=0	
	sp=0
	ge=0
	fa=0
	od=0
	cl=0
	ph=0
	su=0
	for i in range(len(taxLines)):
    		for j in range(i + 1, len(taxLines)):

  		        if taxLines[i].split("\t")[0] != taxLines[j].split("\t")[0]:
				if taxLines[i].split("\t")[1] == taxLines[j].split("\t")[1]:
					sp=sp+1

				elif taxLines[i].split("\t")[1].split(";")[5] == taxLines[j].split("\t")[1].split(";")[5]:
                        		ge=ge+1
				elif taxLines[i].split("\t")[1].split(";")[4] == taxLines[j].split("\t")[1].split(";")[4]:
                        	        fa=fa+1
				elif taxLines[i].split("\t")[1].split(";")[3] == taxLines[j].split("\t")[1].split(";")[3]:
                                        od=od+1
				elif taxLines[i].split("\t")[1].split(";")[2] == taxLines[j].split("\t")[1].split(";")[2]:
                                        cl=cl+1
		        	elif taxLines[i].split("\t")[1].split(";")[1] == taxLines[j].split("\t")[1].split(";")[1]:
                                        ph=ph+1

		        	else:
                                        su=su+1

	print "species: " + str(sp)
	print "genus: "  + str(ge)
	print "family: " + str(fa) 
	print "order: "	+ str(od) 
	print "class: "	+ str(cl) 
	print "phylum: " + str(ph)
	print "sobra: " + str(su)	



# read distances File and write separate files for intra and inter species distances


        disLines = open(disFile, "r").read().splitlines()
        pairsFile = disFile.replace(".dist", "_pairs.dist")

        print "Creating distances file " + pairsFile + "\n"
        paDist = open(pairsFile, "w")

        pairs=[]
        for line in disLines:
                if line.split("\t")[0] != line.split("\t")[1]:
                        pairs.append(''.join((line.split("\t")[0],line.split("\t")[1])))
                        if ''.join((line.split("\t")[1],line.split("\t")[0])) not in pairs:
                       		if accTax[line.split("\t")[0]] == accTax[line.split("\t")[1]]:
                               		paDist.write("SP\t" + accTax[line.split("\t")[0]] + "\t" + line.split("\t")[0] + "\t" + accTax[line.split("\t")[1]] + "\t" + line.split("\t")[1] + "\t" + line.split("\t")[2] + "\n")
                       		elif accTax[line.split("\t")[0]].split(";")[5] == accTax[line.split("\t")[1]].split(";")[5]:
                               		paDist.write("GE\t" + accTax[line.split("\t")[0]] + "\t" + line.split("\t")[0] + "\t" + accTax[line.split("\t")[1]] + "\t" + line.split("\t")[1] + "\t" + line.split("\t")[2] + "\n")
                        	elif accTax[line.split("\t")[0]].split(";")[4] == accTax[line.split("\t")[1]].split(";")[4]:
                                	paDist.write("FA\t" + accTax[line.split("\t")[0]] + "\t" + line.split("\t")[0] + "\t" + accTax[line.split("\t")[1]] + "\t" + line.split("\t")[1] + "\t" + line.split("\t")[2] + "\n")
                        	elif accTax[line.split("\t")[0]].split(";")[3] == accTax[line.split("\t")[1]].split(";")[3]:
                                	paDist.write("OR\t" + accTax[line.split("\t")[0]] + "\t" + line.split("\t")[0] + "\t" + accTax[line.split("\t")[1]] + "\t" + line.split("\t")[1] + "\t" + line.split("\t")[2] + "\n")
                        	elif accTax[line.split("\t")[0]].split(";")[2] == accTax[line.split("\t")[1]].split(";")[2]:
                                	paDist.write("CL\t" + accTax[line.split("\t")[0]] + "\t" + line.split("\t")[0] + "\t" + accTax[line.split("\t")[1]] + "\t" + line.split("\t")[1] + "\t" + line.split("\t")[2] + "\n")
                        	elif accTax[line.split("\t")[0]].split(";")[1] == accTax[line.split("\t")[1]].split(";")[1]:
                                	paDist.write("PH\t" + accTax[line.split("\t")[0]] + "\t" + line.split("\t")[0] + "\t" + accTax[line.split("\t")[1]] + "\t" + line.split("\t")[1] + "\t" + line.split("\t")[2] + "\n")
                        	else:
                               		paDist.write("SU\t" + accTax[line.split("\t")[0]] + "\t" + line.split("\t")[0] + "\t" + accTax[line.split("\t")[1]] + "\t" + line.split("\t")[1] + "\t" + line.split("\t")[2] + "\n")

        paDist.close()


main()
