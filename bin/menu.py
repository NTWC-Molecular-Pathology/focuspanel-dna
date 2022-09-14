import os

ans=True
while ans:
    print ("""
    1. Run Somatic Pipeline
    2. Primer Mapping and SNP Check
    3. ...
    4. Exit/Quit
    """)
    ans=input("What would you like to do?" )
    if ans=="1": 
      os.system("ls ~/lungpanel/input | cut -f 1")
      sample = input ("Please input sample name: ")
      rundate = input ("Please input library run date: ")
      os.system(f'nextflow run ~/lungpanel/nextflow/new_somatic/main.nf --sample {sample} --rundate {rundate}')
    elif ans=="2":
      os.system('conda activate primercheck')
      os.chdir("/home/tmhngs/primercheck/primer-map-snp-check-new_insilico_pcr")
      os.system("python primercheck.py") 
    elif ans=="3":
      print("\n ...") 
    elif ans=="4":
      ans = False
      print("\n Goodbye") 
    elif ans !="":
      print("\n Not Valid Choice Try again") 