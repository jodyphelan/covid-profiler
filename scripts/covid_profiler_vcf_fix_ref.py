import sys
import argparse
import pathogenprofiler as pp


def main(args):
    refseq = list(pp.fasta(args.ref).fa_dict.values())[0]
    for l in sys.stdin:
        if l[0]=="#":
            if l.strip()=="##contig=<ID=1,length=29903>":
                l = l.replace("1","NC_045512.2")
            sys.stdout.write(l)
            continue
        row = l.strip().split()
        if row[0]=="1": row[0] = args.seqname
        ipos = int(row[1])-1
        possible_ref_allele = row[3]
        true_ref_allele = refseq[ipos]
        alts = row[4].split(",")
        alleles = [possible_ref_allele] + alts
        if possible_ref_allele!=true_ref_allele:
            ref_allele_index = alts.index(true_ref_allele)
            row[3] = true_ref_allele
            alts[ref_allele_index] = possible_ref_allele
            row[4] = ",".join(alts)
            genos = "".join(row[9:]).replace(str(ref_allele_index+1),"R").replace("0",str(ref_allele_index+1)).replace("R","0")

            row[9:] = list(genos)
            sys.stdout.write("\t".join(row)+"\n")
        else:
            sys.stdout.write("\t".join(row)+"\n")



parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--ref',help='VCF file',required=True)
parser.add_argument('--seqname',default="NC_045512.2",help='VCF file')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
