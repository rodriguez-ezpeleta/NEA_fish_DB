

#create logfile with all manip that are going to be done
set.logfile(name=Teleo_samples_20210621.log)

#analyse reads (length, ambiguities) and count number of reads per group
summary.seqs(fasta=Teleo_samples_20210621.fasta, processors=12)
count.groups(group=Teleo_samples_20210621.groups)

#dereplicate
unique.seqs(fasta=Teleo_samples_20210621.fasta)
summary.seqs(fasta=Teleo_samples_20210621.unique.fasta, name=Teleo_samples_20210621.names, processors=12)

#align sequences against the 12S rRNA database
align.seqs(fasta=Teleo_samples_20210621.unique.fasta, reference=teleo.align, processors=12, flip=T)
summary.seqs(fasta=Teleo_samples_20210621.unique.align, name=Teleo_samples_20210621.names, processors=12)

#remove sequences not covering the "teleo" region of 12S
screen.seqs(fasta=Teleo_samples_20210621.unique.align, name=Teleo_samples_20210621.names, group=Teleo_samples_20210621.groups, start=6, end=77, minlength=60, maxlength=100, maxambig=0, processors=12)
summary.seqs(fasta=Teleo_samples_20210621.unique.good.align, name=Teleo_samples_20210621.good.names, processors=12)
count.groups(group=Teleo_samples_20210621.good.groups)

# remove columns that contain gap characters
filter.seqs(fasta=Teleo_samples_20210621.unique.good.align, vertical=T, processors=12)
unique.seqs(fasta=Teleo_samples_20210621.unique.good.filter.fasta, name=Teleo_samples_20210621.good.names)
summary.seqs(fasta=Teleo_samples_20210621.unique.good.filter.unique.fasta,name=Teleo_samples_20210621.unique.good.filter.names,processors=12)

#remove chimeras
chimera.uchime(fasta=Teleo_samples_20210621.unique.good.filter.unique.fasta, name=Teleo_samples_20210621.unique.good.filter.names, group=Teleo_samples_20210621.good.groups, processors=12)
remove.seqs(accnos=Teleo_samples_20210621.unique.good.filter.unique.denovo.uchime.accnos, fasta=Teleo_samples_20210621.unique.good.filter.unique.fasta, name=Teleo_samples_20210621.unique.good.filter.names, group=Teleo_samples_20210621.good.groups, dups=T)
summary.seqs(fasta=Teleo_samples_20210621.unique.good.filter.unique.pick.fasta, name=Teleo_samples_20210621.unique.good.filter.pick.names, processors=12)
count.groups(group=Teleo_samples_20210621.good.pick.groups)

#for clarity, rename files
system(cp Teleo_samples_20210621.unique.good.filter.unique.pick.fasta Teleo_samples_20210621_all.fasta)
system(cp Teleo_samples_20210621.unique.good.filter.pick.names Teleo_samples_20210621_all.names)
system(cp Teleo_samples_20210621.good.pick.groups Teleo_samples_20210621_all.groups)

#create count_table with the new files
make.table(name=Teleo_samples_20210621_all.names, group=Teleo_samples_20210621_all.groups, compress=F)

quit

# Clustering and taxonomic assignment
# PHYLOTYPES (assign taxonomy without clustering)
mothur "#classify.seqs(fasta=Teleo_samples_20210621_all.fasta, template=Sequences_Teleo.align, taxonomy=Sequences_Teleo.tax, name=Teleo_samples_20210621_all.names, group=Teleo_samples_20210621_all.groups, method=wang, cutoff=60, processors=12)"
mothur "#phylotype(taxonomy=Teleo_samples_20210621_all.Sequences_Teleo.wang.taxonomy)"
mothur "#make.shared(list=Teleo_samples_20210621_all.Sequences_Teleo.wang.tx.list, count=Teleo_samples_20210621_all.count_table, label=1)"
mothur "#classify.otu(list=Teleo_samples_20210621_all.Sequences_Teleo.wang.tx.list, count=Teleo_samples_20210621_all.count_table, taxonomy=Teleo_samples_20210621_all.Sequences_Teleo.wang.taxonomy, label=1)"
sed -E s/";[a-zA-Z]*_unclassified"/";unclassified"/g Teleo_samples_20210621_all.Sequences_Teleo.wang.tx.1.cons.taxonomy > Teleo_samples_20210621_all.Sequences_Teleo.wang.tx.1.cons_corr.taxonomy

#merge information and create output files
#remove label and phylo_count columns
awk '{$1=$3=""; print $0}' Teleo_samples_20210621_all.Sequences_Teleo.wang.tx.shared > file1
#transpose file
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < file1 > file2
paste file2 Teleo_samples_20210621_all.Sequences_Teleo.wang.tx.1.cons_corr.taxonomy > file3
#check merge is correct
awk '{$4=""; print $0}' file3 > output

rm file*

