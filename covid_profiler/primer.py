import pathogenprofiler as pp
from collections import defaultdict



def parse_fuzznuc_output(output):
    seq = None
    results = defaultdict(list)
    for l in open(output):
        row = l.strip().split()
        if "# Sequence" in l:
            seq = row[2]
            continue
        elif l[0]=="#": continue
        if len(row)<5: continue
        if row[0]=="Start": continue
        results[seq].append({
            "seqname":seq,
            "start":int(row[0]),
            "end":int(row[1]),
            "strand":row[2],
            "pattern":row[3],
            "mismatches":0 if row[4]=="." else int(row[4]),
            "sequence":row[5]
        })
    return results

def run_fuzznuc(seqs,pattern,pmismatch=0):
    tmpfile = pp.get_random_file()
    pp.run_cmd("fuzznuc -sequence %s -pattern %s -outfile %s -complement -pmismatch %s" % (seqs,pattern,tmpfile,pmismatch))
    result = parse_fuzznuc_output(tmpfile)
    # pp.rm_files([tmpfile])
    return result


def find_amplicon(fseqs,rseqs,pseqs):
    if len(fseqs)==1 and len(rseqs)==1 and len(pseqs)==1:
        return {
            "id":fseqs[0]["seqname"],
            "forward_primer_seq":fseqs[0]["sequence"],
            "forward_primer_mismatches":fseqs[0]["mismatches"],
            "probe_sequence":pseqs[0]["sequence"],
            "probe_mismatches":pseqs[0]["mismatches"],
            "reverse_primer_seq":rseqs[0]["sequence"],
            "reverse_primer_mismatches":rseqs[0]["mismatches"],
        }
    else:
        return {
        "id":"NA",
        "forward_primer_seq":"NA",
        "forward_primer_mismatches":"NA",
        "probe_sequence":"NA",
        "probe_mismatches":"NA",
        "reverse_primer_seq":"NA",
        "reverse_primer_mismatches":"NA"
    }
