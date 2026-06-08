# rule merge_assembly_clusters:
#     """
#     Merge MetaFlye contigs + MMseqs representatives.
#     - MetaFlye IDs  : >assembly_<originalID>
#     - MMseqs IDs    : >cluster_<nbReadsInCluster>_<representativeID>
#     """
#     input:
#         assembly = "results/metaflye/{sample}/assembly.fasta",
#         clusters = "results/mmseqs/{sample}/clusters_rep_seq.fasta",
#         summary  = "results/mmseqs/{sample}/clusters_summary.tsv",


def openfile(file):
    """open a file 

    :param file: name of a file
    :type file: string
    :return: list of all the lines of the files
    :rtype: list
    """
    with open(file,'r') as myfile:
        return list(myfile.readlines())
    
def parse_summaryf(summaryf):
    """take summary files and return a dic like dic[ID]=nb reads in cluster
    param summaryf: summary data
    type: list
    return: dic of reads ID and number of reads in the cluster
    rtype: dic
    """
    dic={}
    for l in summaryf:
        tab_ligne=l.split('\t')
        dic[tab_ligne[0]]=tab_ligne[1]
    return dic

def parse_assembly(assemblyf):
    """take assembly files and add labbeled reads ID
    param assemblyf: files data
    type: list
    return: fasta file
    rtype: string
    """
    rep=''
    if len(assemblyf)==0:
        return rep
    else:
        for l in assemblyf:
            if l[0]=='>':
                rep=rep+'>assembly_'+l[1:]
            else:
                rep=rep+l
    return rep 

def parse_mmseq(clusterf,summaryf):
    """take mmseq files and add labbeled reads ID
    param clusterf: files data
    type: list
    param clusterf: summary data
    type: list
    return: fasta file
    rtype: string
    """
    rep=''
    if len(clusterf)==0:
        return rep
    else:
        sum_dic=parse_summaryf(summaryf)
        for l in clusterf:
            if l[0]=='>':
                nb=sum_dic[l.split('\t')[0][1:].strip('\n').strip(' ')]
                rep=rep+'>cluster_'+str(nb)+'_'+l[1:]
            else:
                rep=rep+l
    return rep 

    
def mergeassemblycluster():
    assemblyf=openfile(snakemake.input[0])
    clusterf=openfile(snakemake.input[1])
    summaryf=openfile(snakemake.input[2])
    # assemblyf=openfile("results/metaflye/C26_merged/assembly.fasta")
    # clusterf=openfile("results/mmseqs/C26_merged/clusters_rep_seq.fasta")
    # summaryf=openfile("results/mmseqs/C26_merged/clusters_summary.tsv")

    assembly_res=parse_assembly(assemblyf)
    mmseq_res=parse_mmseq(clusterf,summaryf)
    final_fasta=assembly_res+mmseq_res
    with open(snakemake.output[0],'w') as myfile:
        myfile.write(final_fasta)



mergeassemblycluster()