

#create logfile with all manip that are going to be done
set.logfile(name=MiFish_samples_20210621.log)

#analyse reads (length, ambiguities) and count number of reads per group
summary.seqs(fasta=MiFish_samples_20210621.fasta, processors=12)
count.groups(group=MiFish_samples_20210621.groups)

#dereplicate
unique.seqs(fasta=MiFish_samples_20210621.fasta)
summary.seqs(fasta=MiFish_samples_20210621.unique.fasta, name=MiFish_samples_20210621.names, processors=12)

#align sequences against the 12S rRNA database
align.seqs(fasta=MiFish_samples_20210621.unique.fasta, reference=MiFish.align, processors=12, flip=T)
summary.seqs(fasta=MiFish_samples_20210621.unique.align, name=MiFish_samples_20210621.names, processors=12)

#remove sequences not covering the "MiFish" region of 12S
screen.seqs(fasta=MiFish_samples_20210621.unique.align, name=MiFish_samples_20210621.names, group=MiFish_samples_20210621.groups, minlength=140, maxlength=200, maxambig=0, processors=12)
summary.seqs(fasta=MiFish_samples_20210621.unique.good.align, name=MiFish_samples_20210621.good.names, processors=12)
count.groups(group=MiFish_samples_20210621.good.groups)

# remove columns that contain gap characters
filter.seqs(fasta=MiFish_samples_20210621.unique.good.align, vertical=T, processors=12)
unique.seqs(fasta=MiFish_samples_20210621.unique.good.filter.fasta, name=MiFish_samples_20210621.good.names)
summary.seqs(fasta=MiFish_samples_20210621.unique.good.filter.unique.fasta,name=MiFish_samples_20210621.unique.good.filter.names,processors=12)

#remove chimeras
chimera.uchime(fasta=MiFish_samples_20210621.unique.good.filter.unique.fasta, name=MiFish_samples_20210621.unique.good.filter.names, group=MiFish_samples_20210621.good.groups, processors=12)
remove.seqs(accnos=MiFish_samples_20210621.unique.good.filter.unique.denovo.uchime.accnos, fasta=MiFish_samples_20210621.unique.good.filter.unique.fasta, name=MiFish_samples_20210621.unique.good.filter.names, group=MiFish_samples_20210621.good.groups, dups=T)
summary.seqs(fasta=MiFish_samples_20210621.unique.good.filter.unique.pick.fasta, name=MiFish_samples_20210621.unique.good.filter.pick.names, processors=12)
count.groups(group=MiFish_samples_20210621.good.pick.groups)

#for clarity, rename files
system(cp MiFish_samples_20210621.unique.good.filter.unique.pick.fasta MiFish_samples_20210621_all.fasta)
system(cp MiFish_samples_20210621.unique.good.filter.pick.names MiFish_samples_20210621_all.names)
system(cp MiFish_samples_20210621.good.pick.groups MiFish_samples_20210621_all.groups)

#create count_table with the new files
make.table(name=MiFish_samples_20210621_all.names, group=MiFish_samples_20210621_all.groups, compress=F)

quit()

# Clustering and taxonomic assignment
# PHYLOTYPES (assign taxonomy without clustering)

mothur
classify.seqs(fasta=MiFish_samples_20210621_all.fasta, template=Sequences_MiFish.align, taxonomy=Sequences_MiFish.tax, name=MiFish_samples_20210621_all.names, group=MiFish_samples_20210621_all.groups, method=wang, cutoff=60, processors=12)
phylotype(taxonomy=MiFish_samples_20210621_all.Sequences_MiFish.wang.taxonomy)
make.shared(list=MiFish_samples_20210621_all.Sequences_MiFish.wang.tx.list, count=MiFish_samples_20210621_all.count_table, label=1)
classify.otu(list=MiFish_samples_20210621_all.Sequences_MiFish.wang.tx.list, count=MiFish_samples_20210621_all.count_table, taxonomy=MiFish_samples_20210621_all.Sequences_MiFish.wang.taxonomy, label=1)
quit()
sed -E s/";[a-zA-Z]*_unclassified"/";unclassified"/g MiFish_samples_20210621_all.Sequences_MiFish.wang.tx.1.cons.taxonomy > MiFish_samples_20210621_all.Sequences_MiFish.wang.tx.1.cons_corr.taxonomy

#merge information and create output files
#remove label and phylo_count columns
awk '{$1=$3=""; print $0}' MiFish_samples_20210621_all.Sequences_MiFish.wang.tx.shared > file1
#transpose file
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < file1 > file2
paste MiFish_samples_20210621_all.Sequences_MiFish.wang.tx.1.cons_corr.taxonomy file2 > file3

#check merge is correct
awk '{$4=""; print $0}' file3 > output #OTU table

rm file*
