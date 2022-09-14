import yaml
import sys
import gzip

cov_file = sys.argv[1]

threshold_file = sys.argv[2]

sample_name = cov_file.split(".")[0]

data = dict()


with gzip.open (cov_file, 'rt') as bed:
    for line in bed:
        fields = line.strip().split('\t')
        chr = fields[0]
        start = fields[1]
        end = fields[2]
        pos = "%s:%s-%s"%(chr,start,end)
        amplicon = fields[3]
        gene = amplicon.split("_")[1]
        coverage = float(fields[4])
        #if amplicon.split("_")[1] == "EGFR":
        data[amplicon] = { "Gene": gene, "Chrom": chr, "Start": start, "End": end, "Coverage": coverage }


with gzip.open (threshold_file, 'rt') as bed:
    for line in bed:
        if not line.startswith("#"):
            fields = line.strip().split('\t')
            chr = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            amplicon = fields[3]
            pos = "%s:%s-%s"%(chr,start,end)
            cov500 = float(fields[4])
            pct_500 = (cov500/(end-start))
            data[amplicon]["500x"] = f"{pct_500:.0%}"


header = {
    "id": "cov_amplicon",
    "section_name": "Amplicon Coverage",
    "description": "This plot shows the coverage in each amplicon.",
    "plot_type": "table",
    "data": data,
    'scale': False,
    'col1_header': 'Amplicon',
    'format': '{:,.0f}',
}


with open (f"{sample_name}_mqc.yml", "w") as f:
    yaml.dump(header,f,default_flow_style=False)
