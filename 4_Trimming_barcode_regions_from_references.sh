## Align reference sequences to make reference alignment
mafft --localpair --maxiterate 1000 complete12S_ref_seq.fasta > complete12S_ref_seq.align 
mafft --localpair --maxiterate 1000 completeCOI_ref_seq.fasta > completeCOI_ref_seq.align

## Align available sequences with reference alignments of previous step
mothur "#align.seqs(fasta=Sequences_12S.fasta, reference=complete12S_ref_seq.align, processors=12, flip=T)"
mothur "#align.seqs(fasta=Sequences_COI.fasta, reference=completeCOI_ref_seq.align, processors=12, flip=T)"

# Middle-STEP: remove punctuation marks (dots and dashes) from output sequences
# 12S
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < Sequences_12S.align > Sequences_12S.fasta_lin.align
awk '{gsub(/\.|\-/,"",$2)}1' Sequences_12S.fasta_lin.align > Sequences_12S_trim_temp.fasta
sed -e 's/  */\n/g' Sequences_12S_trim_temp.fasta > Sequences_12S_trimmed.fasta
# COI
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < Sequences_COI.align > Sequences_COI.fasta_lin.align
awk '{gsub(/\.|\-/,"",$2)}1' Sequences_COI.fasta_lin.align > Sequences_COI_trim_temp.fasta
sed -e 's/  */\n/g' Sequences_COI_trim_temp.fasta > Sequences_COI_trimmed.fasta

## Get barcode sequences using primers
# Teleo
cutadapt --discard-untrimmed -g ACACCGCCCGTCACTCT...CATGGTAAGTGTACCGGAAG -o Teleo_12S_ref_seqs.fasta Sequences_12S_trimmed.fasta
mothur "#screen.seqs(fasta=Teleo_12S_ref_seqs.fasta, minlength=50, maxlength=80, processors=12)"
# MiFish
cutadapt --discard-untrimmed -g GTCGGTAAAACTCGTGCCAGC...CAAACTGGGATTAGATACCCCACTATG -o MiFish_12S_ref_seqs.fasta Sequences_12S_trimmed.fasta
mothur "#screen.seqs(fasta=MiFish_12S_ref_seqs.fasta, minlength=150, maxlength=190, processors=12)"
# Leray *using R2 primer
cutadapt --discard-untrimmed -g GGWACWGGWTGAACWGTWTAYCCYCC...TGRTTYTTYGGTCAYCCTGARGTTTA -o Leray_COI_ref_seqs.fasta Sequences_COI_trimmed.fasta
mothur "#screen.seqs(fasta=Leray_COI_ref_seqs.fasta, minlength=300, maxlength=320, processors=12)"
# Folmer
cutadapt --discard-untrimmed -g GGTCAACAAATCATAAAGATATTGG...TGATTTTTTGGTCACCCTGAAGTTTA -o Folmer_COI_ref_seqs.fasta Sequences_COI_trimmed.fasta
mothur "#screen.seqs(fasta=Folmer_COI_ref_seqs.fasta, minlength=640, maxlength=670, processors=12)"

## Align barcode sequences to create barcode reference aligments
## These files will be used in next step as reference sequence file
mafft --localpair --maxiterate 1000 Teleo_12S_ref_seqs_good.fasta > Teleo_12S_ref_seqs.align 
mafft --localpair --maxiterate 1000 MiFish_12S_ref_seqs_good.fasta > MiFish_12S_ref_seqs.align
mafft --localpair --maxiterate 1000 Leray_COI_ref_seqs_good.fasta > Leray_COI_ref_seqs.align 
mafft --localpair --maxiterate 1000 Folmer_COI_ref_seqs_good.fasta > Folmer_COI_ref_seqs.align

## Align available sequences with barcode reference alignments and get all barcode sequences
cp Sequences_12S_trimmed.fasta Sequences_Teleo.fasta
cp Sequences_12S_trimmed.fasta Sequences_MiFish.fasta
cp Sequences_COI_trimmed.fasta Sequences_Leray.fasta
cp Sequences_COI_trimmed.fasta Sequences_Folmer.fasta
mothur "#align.seqs(fasta= Sequences_Teleo.fasta, reference= Teleo_12S_ref_seqs.align, processors=12, flip=T)"
mothur "#align.seqs(fasta= Sequences_MiFish.fasta, reference= MiFish_12S_ref_seqs.align, processors=12, flip=T)"
mothur "#align.seqs(fasta= Sequences_Leray.fasta, reference= Leray_COI_ref_seqs.align, processors=12, flip=T)"
mothur "#align.seqs(fasta= Sequences_Folmer.fasta, reference= Folmer_COI_ref_seqs.align, processors=12, flip=T)"

## Keep sequences with > 90% length 
mothur "#screen.seqs(fasta=Sequences_Folmer.align, minlength=592, maxambig=0, processors=12)" # 592=90% de 658 
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < Sequences_Folmer.good.align > Sequences_Folmer.good_lin.align
awk -F "\t" '{gsub(/\.|\-/,"",$2);print $1"\t"$2}' Sequences_Folmer.good_lin.align > Sequences_Folmer_trim_temp
sed -e 's/\t/\n/g' Sequences_Folmer_trim_temp > Sequences_Folmer.fasta

mothur "#screen.seqs(fasta=Sequences_Leray.align, minlength=281, maxambig=0, processors=12)" # 281 = 90% de 313
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < Sequences_Leray.good.align > Sequences_Leray.good_lin.align
awk -F "\t" '{gsub(/\.|\-/,"",$2);print $1"\t"$2}' Sequences_Leray.good_lin.align > Sequences_Leray_trim_temp
sed -e 's/\t/\n/g' Sequences_Leray_trim_temp > Sequences_Leray.fasta

mothur "#screen.seqs(fasta=Sequences_MiFish.align, minlength=146, maxambig=0, processors=12)" # 146 = 90% de 163-185 
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < Sequences_MiFish.good.align > Sequences_MiFish.good_lin.align
awk -F "\t" '{gsub(/\.|\-/,"",$2);print $1"\t"$2}' Sequences_MiFish.good_lin.align > Sequences_MiFish_trim_temp
sed -e 's/\t/\n/g' Sequences_MiFish_trim_temp > Sequences_MiFish.fasta

mothur "#screen.seqs(fasta=Sequences_Teleo.align, minlength=56, maxambig=0, processors=12)" # 56 = 90% de 63 
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < Sequences_Teleo.good.align > Sequences_Teleo.good_lin.align
awk -F "\t" '{gsub(/\.|\-/,"",$2);print $1"\t"$2}' Sequences_Teleo.good_lin.align > Sequences_Teleo_trim_temp
sed -e 's/\t/\n/g' Sequences_Teleo_trim_temp > Sequences_Teleo.fasta




