#! /usr/bin/env python2.7

'''
This script takes FlyBase ids from the supplementary material of Obbard 2009, together
with their associated Gene Ontology annotations, and for each one determines a list of
orthologous H. melpomene loci.

'''

import csv
import re
import urllib

### DATA FILES
baseDir = "/Users/joe/Arch/JigginsRotation/HelPopGen"
# csv files from Obbard 2009 supplementary material
fAll = "/Users/joe/Arch/JigginsRotation/papers/ObbardEtAl2009_Supplementary_Material/s027_all.csv"
fKenyan = "/Users/joe/Arch/JigginsRotation/papers/ObbardEtAl2009_Supplementary_Material/s027_kenyan.csv"
# orthomcl output file
fOrthomcl = "/Users/joe/Arch/JigginsRotation/HelPopGen/data/all_orthomcl.out.parsed"

def FBtoAnnID(FBgn):
    '''
    Takes a FlyBase id from the supplementary file s027 from Obbard et al 2009, and 
    retreives corresponding Annotation Ids by scraping flybase pages. The
    Annotation Ids match the dataset used for finding orthologues between D 
    melanogaster and H melpomene genes by Will Palmer.
    '''
    try:
        url = "http://flybase.org/reports/" + FBgn + ".html"
        f = urllib.urlopen(url)
        page = f.readlines()
        if f.geturl() == url:
            value = re.findall(
                    r'(?<=Annotation symbol\</th\>\<td\>)[A-Za-z0-9]+',
                    ''.join(page))
            print FBgn, ':  ', value
            return value
        else:
            names = list(set(re.findall(r'(?<=value=")FBgn[0-9]+', ''.join(page))))
            ids = list()
            for name in names:
                ids.extend(FBtoAnnID(name))
            print FBgn, ':  ', ids
            return ids
    except Exception as e:
        f.close()
    else:
        f.close()

def DmeltoHmel(name):
    '''
    Given an Annotation Id (see the FBtoAnnID method), a list of corresponding 
    H melpomene genes is retreived from an orthomcl output file.

    '''
    group = None
    with open(fOrthomcl, 'r') as f:
        for line in f.readlines():
            if name == line.split()[1].split('-')[0]:
                group = line.split()[0]
    if group is None:
        return list()
    HmelGenes = list()
    with open(fOrthomcl, 'r') as f:
        for line in f.readlines():
            if group in line:
                if line.split()[2] == "Hmel":
                    HmelGenes.append(line.split()[1].split('-')[0])
    print HmelGenes
    return HmelGenes

def main():
    with open('s027_all.csv', 'r') as fAll:
        f = csv.reader(fAll)
        colNames = f.next()
        rows = [{k:v for (k,v) in zip(colName, line)} for line in f]
    with open("classes.csv", 'wb') as fResults:
        resultsWriter = csv.writer(fResults)
        headers = ["Hmel_Locus",
                   "Dmel_Annotation_ID",
                   "FlyBase_ID",
                   "Dmel_Locus",
                   "Class",
                   "Cell_Hum",
                   "Ls",
                   "Ln",
                   "Dn",
                   "Ds",
                   "Dmel_Pn",
                   "Dsim_Pn",
                   "Dmel_Ps",
                   "Dsim_Ps",]
        resultsWriter.writerow(headers)
        for row in rows:
            for DmelGene in FBtoAnnID(row["FBgn"]):
                for HmelGene in DmeltoHmel(DmelGene):
                    outRow = [HmelGene,
                              DmelGene,
                              row["FBgn"],
                              row["Locus"],
                              row["Class"],
                              row["Cell_Hum"],
                              row["Ls"],
                              row["Ln"],
                              row["Dn"],
                              row["Ds"],
                              row["Mel_Pn"],
                              row["Sim_Pn"],
                              row["Mel_Ps"],
                              row["Sim_Ps"],]
                    resultsWriter.writerow(outRow)

if __name__ == "__main__":
    main()

