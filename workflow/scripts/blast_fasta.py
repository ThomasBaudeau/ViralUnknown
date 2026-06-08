
def openfile(file):
    """open a file 

    :param file: name of a file
    :type file: string
    :return: list of all the lines of the files
    :rtype: list
    """
    with open(file,'r') as myfile:
        return list(myfile.readlines())
    
def parse_candidatef(candidatef):
    """take multifasta files and return a dic like dic[ID]=fasta
    param summaryf: multifasta data
    type: list
    return: dic of reads ID /fasta
    rtype: dic
    """
    dic={}
    id=None
    for l in candidatef:
        if l[0]=='>':
            if id:
                dic[id]=rep
            id=l.split('\t')[0][1:].strip('\n').strip(' ')
            rep=l
        else:
            rep=rep+l
    dic[id]=rep
    return dic

def parse_blastsum(blastsumf):
    """take blast summary files and extract if reads need to be keept
    param assemblyf: files data
    type: list
    return: lists of discarded or saved reads
    rtype: list
    """
    keep=[]
    discard=[]
    for l in blastsumf[1:]:
        tab_line=l.split('\t')
        if tab_line[1]=='excluded':
            discard.append(tab_line[0])
        else:
            keep.append(tab_line[0])
    return keep,discard 



    
def mergeassemblycluster():
    candidatef=openfile(snakemake.input[0])
    blastsumf=openfile(snakemake.input[1])
    # candidatef=openfile("results/blast_input/C26_merged_candidates.fasta")
    # blastsumf=openfile("results/blast/C26_merged_blast_summary.tsv")

    dic_candidat=parse_candidatef(candidatef)
    keep,discard=parse_blastsum(blastsumf)

    # with open('C26_merged_candidates_filtered.fasta','w') as filtered:
    #     for id in keep:
    #         filtered.write(dic_candidat[id])
    # with open('C26_merged_candidates_excluded.fasta','w') as excluded:
    #     for id in discard:
    #         excluded.write(dic_candidat[id])

    print(dic_candidat.keys())
    with open(snakemake.output[0],'w') as filtered:
        for id in keep:
            filtered.write(dic_candidat[id])
    with open(snakemake.output[1],'w') as excluded:
        for id in discard:
            excluded.write(dic_candidat[id])

mergeassemblycluster()