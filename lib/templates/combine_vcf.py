import os
import pandas as pd
import sys

sample_name = sys.argv[1]


cwd = os.path.abspath('')
files = os.listdir(cwd)


df = pd.DataFrame()

file_name = f'{sample_name}.combine_vcf.xlsx'

with pd.ExcelWriter(file_name) as writer:
    for file in files:
        if file.endswith('.txt'):
            if '1000x' in file:
                sheetname = f'mutect2 1000x'
            elif 'mutect' in file:
                sheetname = 'mutect2'
            else:
                sheetname = 'freebayes'
            df = pd.read_csv(file, sep='\t', header= None)
            df.to_excel(writer, sheet_name = sheetname, index = False, header = None )







