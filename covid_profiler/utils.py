import pathogenprofiler as pp
from collections import Counter, defaultdict
import csv
import re
import sys

def get_conf_dict(library_prefix):
    files = {"gff":".gff","ref":".fasta","barcode":".barcode.bed","version":".version.json","proteins":".proteins.csv","msa":".msa.fa","non_coding_bed":".non_coding.bed","meta":".msa.meta.csv"}
    conf = {}
    for key in files:
        sys.stderr.write("Using %s file: %s\n" % (key,library_prefix+files[key]))
        conf[key] = pp.filecheck(library_prefix+files[key])
    return conf

def log(string):
    sys.stderr.write("\n"+"*" * 40+"\n")
    sys.stderr.write(string+"\n")
    sys.stderr.write("*" * 40+"\n\n")

def get_codon_num(mutation):
    re_obj = re.search("([\d]+)[A-Z\*]",mutation)
    if re_obj:
        return int(re_obj.group(1))
    else:
        return None

def get_N_content(seq):
    cnt = Counter(seq)
    return cnt["N"]/len(seq)


def change_codon_number(mut,num):
    re_obj = re.search("(\d+)([A-Z\*]+)>(\d+)([A-Z\*]+)",mut)
    if re_obj:
        return "%s%s>%s%s" % (num,re_obj.group(2),num,re_obj.group(4))
    re_obj = re.search("(\d+)([A-Z\*]+)",mut)
    if re_obj:
        return "%s%s" % (num,re_obj.group(2))
    else:
        import pdb; pdb.set_trace()


def get_variant_data(vcf_file,ref_file,gff_file,protein_file):
    nsp_data = {}
    gene_info = {}
    for row in csv.DictReader(open(protein_file)):
        row["Start"] = int(row["Start"])
        row["End"] = int(row["End"])
        gene_info[row["Gene"]] = {"function":row["Putative function"],"DOI":row["DOI"]}
        if row["Region"]!="nsp": continue
        for i in range(row["Start"],row["End"]+1):
            nsp_data[i] = row

    pp.run_cmd("samtools faidx %s" % ref_file)
    results = defaultdict(list)
    for l in pp.cmd_out("bcftools view %s | bcftools csq -f %s -g %s -p a  | correct_covid_csq.py | bcftools +fill-tags | bcftools query -f '%%POS\\t%%REF\\t%%ALT\\t%%AF\\t%%BCSQ\\n'" % (vcf_file,ref_file,gff_file)):
        # Replace " " with "N" because if the alt allele contains N then
        # translated consequence will have spaces
        row = l.strip().replace(" ","N").split()
        pos,ref,alts_str,af_str,csq_str = row
        alt_af = sum([float(x) for x in af_str.split(",")])
        csqs = csq_str.split(",")
        types = []
        changes = []
        genes = []
        pos = int(pos)
        for i in range(len(csqs)):
            if csqs[i][0]=="@":
                # results[pos].append(results[int(csqs[i][1:])][0])
                pass
            elif csqs[i]==".":
                results[pos].append({"pos":pos, "alts":alts_str, "alt_af":alt_af, "types":"intergenic","changes":"NA","gene":"NA","gene_function":"NA","gene_reference":"NA"})

            else:
                csq = csqs[i].split("|")
                types.append(csq[0].replace("*",""))
                if csq[1]=="orf1ab":
                    codon_pos = get_codon_num(csq[5])
                    if codon_pos in nsp_data:
                        genes.append(nsp_data[codon_pos]["Gene"])
                        codon_pos = codon_pos-nsp_data[codon_pos]["Start"]+1
                        changes.append(change_codon_number(csq[5],codon_pos))
                    else:
                        genes.append("orf1ab")
                        changes.append(csq[5])
                else:
                    changes.append(csq[5] if len(csq)>5 else "")
                    genes.append(csq[1])
                if len(set(types))==1:
                    types = list(set(types))
                results[pos].append({"pos":pos, "alts":alts_str, "alt_af":alt_af, "types":",".join(types), "changes":",".join(changes),"gene":genes[0], "gene_function":gene_info[genes[0]]["function"], "gene_reference":gene_info[genes[0]]["DOI"]})
    final_results = []
    for res in list(results.values()):
        for r in res:
            final_results.append(r)
        # if len(res)==1:
        #     final_results.append(res[0])
        # else:
        #     quit("ERROR! more than one variant for a position")
    return final_results


def vcf2consensus(bam,vcf,ref,id,consensus):
    tmp_bed = pp.get_random_file()
    pp.run_cmd('bedtools genomecov -d -ibam %s | awk \'$3<10\' | awk \'{print $1"\\t"$2"\\t"$2}\' > %s' % (bam,tmp_bed))
    pp.run_cmd("bcftools consensus -f %s -m %s -M N %s | sed 's/^>.*/>%s/' > %s" % (ref,tmp_bed,vcf,id,consensus))
    pp.run_cmd("rm %s" % tmp_bed)
