import glob, tabula
import pandas as pd

pdfs = glob.glob("*.pdf")

for i in pdfs:
    name=i.split('.')[0].replace(" ","")
    print("starting " + name)
    dfs = tabula.read_pdf(i, pages='2')
    df=dfs[0]
    df = df.rename(columns={'Unnamed: 0': 'number', 'Cytogenetic\rLocation': 'cytogenetic_location', 
                            'Start': 'Start_grch37', 'Stop': 'Stop_grch37', 'Gain/\rLoss': 'Type', 
                            'Classification (Initial)': 'Classification_initial', 'Classification (Final)': 'Classification_final'})
    df=df.replace({'\r': ''}, regex=True)
    df.to_csv(name+"_array_cnvs.csv", index=False)
