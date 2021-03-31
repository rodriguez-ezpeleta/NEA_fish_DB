# from a file containing a list of species name, retrieve from Genebank "Species_Name Taxa_Id Gi Acc_Number Sequence_Title Sequence_Length"

# usage python getGenes.py sp_names.txt sp_genes.txt

def main():

    import os
    import sys
    import string
    import re
    import urllib

    eUtilsUrl = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    eSearchUrl = 'esearch.fcgi?'
    eLinkUrl = 'elink.fcgi?'
    eFetchUrl = 'efetch.fcgi?'
    eSummaryUrl = 'esummary.fcgi?'

    taxFile = sys.argv[1]
    genFile = sys.argv[2]
    inFile = open(taxFile, 'r')
    outFile = open(genFile, 'w')

    taxs = inFile.readlines()

    regexpId =  re.compile(r"(<.*Id.*>)"r"(?P<Id>.+)"r"</Id>")
    regexpTax = re.compile(r"<Lineage>"r"(?P<Lineage>.+)"r"</Lineage>")
    regexpSp =  re.compile(r"<ScientificName>"r"(?P<ScientificName>.+)"r"</ScientificName>")
    regexpTitle = re.compile(r"((<Item.*(Title)+.*>))"r"(?P<Title>.+)"r"((<\/Item>))")
    regexpTaxId = re.compile(r"((<Item.*(TaxId)+.*>))"r"(?P<TaxId>.+)"r"((<\/Item>))")
    regexpLength = re.compile(r"((<Item.*(Length)+.*>))"r"(?P<Length>.+)"r"((<\/Item>))")
    regexpAcc = re.compile(r"((<Item.*(AccessionVersion)+.*>))"r"(?P<Acc>.+)"r"((<\/Item>))")
    outFile.write("Species_Name\tTaxa_Id\tGi\tAcc_Number\tSequence_Title\tSequence_Length\n")


    for tax in taxs:

        eSearchQuery = eUtilsUrl + eSearchUrl + 'db=nucleotide&retmax=100000&term=' + tax + "[orgn] AND gene_in_mitochondrion[prop]"
        accessionUidTable = urllib.urlopen(eSearchQuery).readlines()
        for line in accessionUidTable :
            result = regexpId.search(line)
            if result :
                accId = result.group('Id')
                eSummaryQuery = eUtilsUrl + eSummaryUrl + 'db=nucleotide&id=' + accId
                genData = urllib.urlopen(eSummaryQuery).readlines()

                for line in genData :
                    isId = regexpId.search(line)
                    isTax = regexpTaxId.search(line)
                    isTitle = regexpTitle.search(line)
                    isLength = regexpLength.search(line)
                    isAcc = regexpAcc.search(line)
                    if isId:
                        Id = isId.group('Id')
                    if isTax:
                        TaxId = isTax.group('TaxId')
                    if isTitle:
                        Title = isTitle.group('Title')
                    if isLength:
                        Length = isLength.group('Length')
                    if isAcc:
                        Acc = isAcc.group('Acc')

                outFile.write(tax.strip() + "\t" + TaxId  + "\t" + Id  + "\t" + Acc + "\t" + Title + "\t" + Length + "\n")

    inFile.close()
    outFile.close()

main()
