import sys
import pyfastx
import argparse
from tqdm import tqdm
import pathogenprofiler as pp
from collections import defaultdict,Counter
import csv
import json
import math
import covid_profiler
import seqlogo
import uuid
import pandas as pd
import plotly.express as px
from urllib.request import urlopen
import re
import os


def revcom(s):
        """Return reverse complement of a sequence"""
        def complement(s):
                        basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',"-":"-"}
                        letters = list(s)
                        letters = [basecomplement[base] for base in letters]
                        return ''.join(letters)
        return complement(s[::-1])

def blast_primer(ref,primer,nident_fraction=0.7):
    primer_fasta = str(uuid.uuid4())
    blast_output = str(uuid.uuid4())

    with open(primer_fasta,"w") as O:
        O.write(">primer_seq\n%s\n" % (primer))
    ident_cutoff = math.ceil(len(primer)*nident_fraction)
    pp.run_cmd("blastn -task blastn-short -query %s -subject %s -outfmt '6 sstart send sseq length mismatch qstart qend nident' > %s" % (primer_fasta,ref,blast_output))
    if open(blast_output).readline().strip()=="":
        return None
    row = open(blast_output).readline().strip().split()
    start,end,matchseq,matchlen,mismatch,qstart,qend,nident = int(row[0]),int(row[1]),row[2],int(row[3]),int(row[4]),int(row[5]),int(row[6]),int(row[7])
    if nident<ident_cutoff:
        return None
    start = start - (qstart-1)
    end = end + (len(primer) - qend)
    pp.rm_files([primer_fasta,blast_output])
    return {"start":start,"end":end,"matchseq":matchseq,"mismatch":mismatch,"matchlen":matchlen}


acgt = set(["A","C","G","T",""])

iupac = {
    "A": set(["A"]),
    "C": set(["C"]),
    "G": set(["G"]),
    "T": set(["T"]),
    "R": set(["A", "G"]),
    "Y": set(["C", "T"]),
    "S": set(["G", "C"]),
    "W": set(["A", "T"]),
    "K": set(["G", "T"]),
    "M": set(["A", "C"]),
    "B": set(["C", "G", "T"]),
    "D": set(["A", "G", "T"]),
    "H": set(["A", "C", "T"]),
    "V": set(["A", "C", "G"]),
    "N": set(["A", "C", "G", "T"])
}

def extract_region_from_fasta(fasta,start,end):
    (p1,p2) = (start,end) if start<end else (end,start)
    seqs = []
    for entry in tqdm(pyfastx.Fasta(fasta, full_name=True)):
        tmp_seq = entry.seq[p1-1:p2].upper()
        if start>end:
            try:
                tmp_seq = revcom(tmp_seq)
            except:
                print(entry.name)
                print(tmp_seq)
                quit()
        seqs.append({"name":entry.name,"seq":tmp_seq})
    return seqs

def get_mismatches(ref,primer,fasta):
    blast_result = blast_primer(ref,primer)
    print(blast_result)
    seqs = extract_region_from_fasta(fasta,blast_result["start"],blast_result["end"])
    mismatches = []
    for seq in seqs:
         num_snps = sum([1 for x,y in zip(primer,seq["seq"]) if x in acgt and y in iupac and x not in iupac[y]])
         num_missing = sum([1 for x in seq["seq"] if x not in acgt])
         mismatches.append({"id":seq["name"],"num_snps":num_snps,"num_missing":num_missing,"seq":seq["seq"]})
    return mismatches

def get_sequence_logo(seqs):

    count_pm = [
        [float(0) for _ in range(len(seqs[0]))],
        [float(0) for _ in range(len(seqs[0]))],
        [float(0) for _ in range(len(seqs[0]))],
        [float(0) for _ in range(len(seqs[0]))],
    ]

    for i in tqdm(range(len(seqs))):
        for j in range(len(seqs[0])):
            for k,n in enumerate(["A","C","G","T"]):
                if seqs[i][j]==n:
                    count_pm[k][j]+=1


    pm = pd.DataFrame(data={"A":count_pm[0],"C":count_pm[1],"G":count_pm[2],"T":count_pm[3]})
    ppm = pm.apply(axis=1,func=lambda x: [y/sum(x) for y in x],result_type='broadcast')
    print(ppm)
    cpm = seqlogo.CompletePm(ppm = ppm)
    tmp_file = "/tmp/%s.svg" % uuid.uuid4()
    seqlogo.seqlogo(cpm, ic_scale = False, format = 'svg', size = 'medium',filename="%s" %tmp_file)
    svg = open(tmp_file).read()
    os.remove(tmp_file)
    for x in set(re.findall("glyph[0-9]+-[0-9]+",svg)):
        svg = svg.replace(x,str(uuid.uuid4()))

    return svg


def main(args):
    conf = covid_profiler.get_conf_dict(sys.base_prefix+"/share/covidprofiler/%s" % args.db)
    vars(args).update(conf)

    if not args.fasta:
        args.fasta = conf["msa"]
    if not args.meta:
        args.meta = conf["meta"]
    if not args.ref:
        args.ref = conf["ref"]
    args.out = "%s/%s" % (args.dir,args.out)

    fp_snps = get_mismatches(args.ref,args.fp,args.fasta)
    probe_snps = get_mismatches(args.ref,args.probe,args.fasta)
    rp_snps = get_mismatches(args.ref,args.rp,args.fasta)

    isolates = [x["id"] for x in fp_snps]
    results = {}
    for res in fp_snps:
        results[res["id"]] = {"id":res["id"]}
        for x in res:
            results[res["id"]]["fp_"+x] = res[x]

    for res in probe_snps:
        for x in res:
            if x=="id": continue
            results[res["id"]]["probe_"+x] = res[x]

    for res in rp_snps:
        for x in res:
            if x=="id": continue
            results[res["id"]]["rp_"+x] = res[x]

    if args.meta:
        isolates_in_meta = set()
        for row in csv.DictReader(open(args.meta)):
            if row["id"] in results:
                isolates_in_meta.add(row["id"])
                for x in row:
                    results[row["id"]][x] = row[x]

    results = [r for r in results.values() if r["id"] in isolates_in_meta] if args.meta else list(results.values())
    # print(results[0])
    if args.write_csv:
        with open(args.out+".primer_results.csv","w") as O:
            writer = csv.DictWriter(O,fieldnames=list(results[0]))
            writer.writeheader()
            writer.writerows(results)


    json_results = {
        "fp_global_mismatch":len([x for x in results if x["fp_num_snps"]>0])/len(results),
        "rp_global_mismatch":len([x for x in results if x["rp_num_snps"]>0])/len(results),
        "probe_global_mismatch":len([x for x in results if x["probe_num_snps"]>0])/len(results)
    }




    with urlopen('https://raw.githubusercontent.com/jodyphelan/tbdr/main/tbdr/static/custom.geo.json') as response:
        geojson = json.load(response)

    csv_df = pd.DataFrame(data={
        "fp_num_snps": [x["fp_num_snps"] for x in results],
        "rp_num_snps": [x["rp_num_snps"] for x in results],
        "probe_num_snps": [x["probe_num_snps"] for x in results],
        "continent": [x["continent"] for x in results],
        "date": [x["date"] for x in results],
        "iso_a3": [x["iso_a3"] for x in results],
    })
    # csv_df = pd.read_csv(args.out+".primer_results.csv")
    csv_df["date"] = pd.to_datetime(csv_df["date"])
    csv_df.set_index('date',inplace=True)

    freq_continent_time = defaultdict(dict)

    for continent in csv_df.continent.unique():
        freq_continent_time["fp_num_snps"][continent] = csv_df[csv_df["continent"]==continent]["fp_num_snps"].resample("1M").apply(lambda x: sum([1 for d in x if d>0])/len(x)*100 if len(x)>0 else None)
        freq_continent_time["rp_num_snps"][continent] = csv_df[csv_df["continent"]==continent]["rp_num_snps"].resample("1M").apply(lambda x: sum([1 for d in x if d>0])/len(x)*100 if len(x)>0 else None)
        freq_continent_time["probe_num_snps"][continent] = csv_df[csv_df["continent"]==continent]["probe_num_snps"].resample("1M").apply(lambda x: sum([1 for d in x if d>0])/len(x)*100 if len(x)>0 else None)

    freq_continent_time["fp_num_snps"]["Total"] = csv_df["fp_num_snps"].resample("1M").apply(lambda x: sum([1 for d in x if d>0])/len(x)*100 if len(x)>0 else None)
    freq_continent_time["rp_num_snps"]["Total"] = csv_df["rp_num_snps"].resample("1M").apply(lambda x: sum([1 for d in x if d>0])/len(x)*100 if len(x)>0 else None)
    freq_continent_time["probe_num_snps"]["Total"] = csv_df["probe_num_snps"].resample("1M").apply(lambda x: sum([1 for d in x if d>0])/len(x)*100 if len(x)>0 else None)

    freq_continent_time["fp_num_snps"] = pd.DataFrame(freq_continent_time["fp_num_snps"])
    freq_continent_time["rp_num_snps"] = pd.DataFrame(freq_continent_time["rp_num_snps"])
    freq_continent_time["probe_num_snps"] = pd.DataFrame(freq_continent_time["probe_num_snps"])

    fig_fp_time = px.line(freq_continent_time["fp_num_snps"], y=freq_continent_time["fp_num_snps"].columns)
    fig_fp_time.update_xaxes(dtick="M1",tickformat="%b\n%Y")


    df = csv_df.pivot_table(values=["fp_num_snps","probe_num_snps","rp_num_snps"], index="iso_a3", aggfunc=lambda x: sum([1 for d in x if d>0])/len(x)*100)
    df["iso_a3"]  = df.index

    fig_fp_map = px.choropleth(
                            df, geojson=geojson, locations='iso_a3', color='fp_num_snps',
                            color_continuous_scale="Viridis",featureidkey='properties.iso_a3',
                            range_color=(0, 100),
                            scope="world",
                            labels={'fp_num_snps':'% samples with SNP'}
                       )
    fig_fp_time = px.line(freq_continent_time["fp_num_snps"], y=freq_continent_time["fp_num_snps"].columns,template="simple_white")
    fig_fp_time.update_xaxes(dtick="M1",tickformat="%b\n%Y")

    json_results["fp_time"] = fig_fp_time.to_html(include_plotlyjs=False,full_html=False)
    json_results["fp_map"] = fig_fp_map.to_html(include_plotlyjs=False,full_html=False)

    fig_rp_map = px.choropleth(
                            df, geojson=geojson, locations='iso_a3', color='rp_num_snps',
                            color_continuous_scale="Viridis",featureidkey='properties.iso_a3',
                            range_color=(0, 100),
                            scope="world",
                            labels={'fp_num_snps':'% samples with SNP'},
                        )
    fig_rp_time = px.line(freq_continent_time["rp_num_snps"], y=freq_continent_time["rp_num_snps"].columns,template="simple_white")
    fig_rp_time.update_xaxes(dtick="M1",tickformat="%b\n%Y")

    json_results["rp_time"] = fig_rp_time.to_html(include_plotlyjs=False,full_html=False)
    json_results["rp_map"] = fig_rp_map.to_html(include_plotlyjs=False,full_html=False)

    fig_probe_map = px.choropleth(
                            df, geojson=geojson, locations='iso_a3', color='probe_num_snps',
                            color_continuous_scale="Viridis",featureidkey='properties.iso_a3',
                            range_color=(0, 100),
                            scope="world",
                            labels={'fp_num_snps':'% samples with SNP'},
                        )
    fig_probe_time = px.line(freq_continent_time["probe_num_snps"], y=freq_continent_time["probe_num_snps"].columns,template="simple_white")
    fig_probe_time.update_xaxes(dtick="M1",tickformat="%b\n%Y")

    json_results["probe_time"] = fig_probe_time.to_html(include_plotlyjs=False,full_html=False)
    json_results["probe_map"] = fig_probe_map.to_html(include_plotlyjs=False,full_html=False)


    json_results["fp_seqlogo"] = get_sequence_logo([x["fp_seq"] for x in results])
    json_results["rp_seqlogo"] = get_sequence_logo([x["rp_seq"] for x in results])
    json_results["probe_seqlogo"] = get_sequence_logo([x["probe_seq"] for x in results])
    # import pdb; pdb.set_trace()
    if args.write_json:
        json.dump(json_results,open(args.out+".plots.json","w"))
    else:
        with open(args.out+".plots.html","w") as O:
            O.write("""
                <html>
                    <head>
                        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
                    </head>
                    <body>
                    <h1>Forward Primer</h1>
                    <h2>Sequence Logo</h2>
                    %(fp_seqlogo)s
                    <h2>Map</h2>
                    %(fp_map)s
                    <h2>Time<h2>
                    %(fp_time)s

                    <h1>Reverse Primer</h1>
                    <h2>Sequence Logo</h2>
                    %(rp_seqlogo)s
                    <h2>Map</h2>
                    %(rp_map)s
                    <h2>Time<h2>
                    %(rp_time)s

                    <h1>Probe</h1>
                    <h2>Sequence Logo</h2>
                    %(probe_seqlogo)s
                    <h2>Map</h2>
                    %(probe_map)s
                    <h2>Time<h2>
                    %(probe_time)s

                    </body>
                </html>
            """ % json_results)



parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--fasta',type=str,help='File with samples')
parser.add_argument('--ref',type=str,help='File with samples')
parser.add_argument('--fp',help='File with samples',required=True)
parser.add_argument('--probe',help='File with samples',required=True)
parser.add_argument('--rp',help='File with samples',required=True)
parser.add_argument('--meta',help='File with samples')
parser.add_argument('--out',help='File with samples',required=True)
parser.add_argument('--plots',action="store_true",help='File with samples')
parser.add_argument('--write-csv',action="store_true",help='File with samples')
parser.add_argument('--write-json',action="store_true",help='File with samples')
parser.add_argument('--db',default="cvdb",help='File with samples')
parser.add_argument('--dir',default=".",help='File with samples')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
