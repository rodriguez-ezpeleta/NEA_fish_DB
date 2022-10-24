## ADDING NEW SEQUENCES TO REFERENCE DATABASE ##

#Directory: /share/projects/OCEAN_eDNA/Cristina/Genes/6_New_sequences_to_refDB/

# Required files: 
	# New sequences file (fasta)
	# Barcode reference alignments (created in step 4_trimming)
	cp Teleo_12S_ref_seqs.align ../6_New_sequences_to_refDB/
	cp MiFish_12S_ref_seqs.align ../6_New_sequences_to_refDB/
	
	#0. Remove marker sequences from fasta with BioEdit
	#0. Rename the fasta with new sequeces to match the reference db (.R)
	
	#1. Align available sequences with reference alignments 
mothur "#align.seqs(fasta=New_sequences_12S_renamed.fasta, reference=Teleo_12S_ref_seqs.align, processors=12, flip=T)"
mv New_sequences_12S_renamed.align New_sequences_teleo.align 
mothur "#align.seqs(fasta=New_sequences_12S_renamed.fasta, reference=MiFish_12S_ref_seqs.align, processors=12, flip=T)"
mv New_sequences_12S_renamed.align New_sequences_mifish.align
	#2. Turn align into fasta
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < New_sequences_teleo.align > New_sequences_teleo_lin.align
awk '{gsub(/\.|\-/,"",$2)}1' New_sequences_teleo_lin.align > New_sequences_teleo_trim_temp.fasta
sed -e 's/  */\n/g' New_sequences_teleo_trim_temp.fasta > New_sequences_teleo_trimmed.fasta
	
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < New_sequences_mifish.align > New_sequences_mifish_lin.align
awk '{gsub(/\.|\-/,"",$2)}1' New_sequences_mifish_lin.align > New_sequences_mifish_trim_temp.fasta
sed -e 's/  */\n/g' New_sequences_mifish_trim_temp.fasta > New_sequences_mifish_trimmed.fasta

	#3. Add aligned sequences into curated reference databases
cat New_sequences_teleo_trimmed.fasta Sequences_Teleo_curated.fasta > Sequences_Teleo_curated_newseqs.fasta
cat New_sequences_mifish_trimmed.fasta Sequences_MiFish_curated.fasta > Sequences_MiFish_curated_newseqs.fasta
	#4. Prepare tax file (.R)
	#5. Do taxa assignment
cp Sequences_Teleo_curated_newseqs.fasta ../Teleo_MiFish/Teleo/Tax_assign/DB_curated/
cp Sequences_MiFish_curated_newseqs.fasta ../Teleo_MiFish/Mifish/Tax_assign/DB_curated/

cd ../Teleo_MiFish/Teleo/Tax_assign/DB_curated/
cd ../Teleo_MiFish/Mifish/Tax_assign/DB_curated/