import pandas as pd
import pyranges as pr
import glob

def alamut_to_pyranges(alamut_df):
    ''' Convert fakealamut file into a pyranges object and
        remove duplicate CNVs (b/c of multiple transcripts etc) '''
    # re-format df to be in a bed-like format for pyranges
    pyranges_df = alamut_df[['chrom','gDNAstart','gDNAend','gene','type']]
    pyranges_df.rename({'chrom': 'Chromosome', 'gDNAstart': 'Start', 'gDNAend': 'End'}, axis=1, inplace=True)
    # convert df to pyranges object
    pyranges = pr.PyRanges(df=pyranges_df)
    # remove duplicate CNVs (these will be because of multiple transcripts). slack -1 i think makes it inclusive of the final position?
    print("Removing duplicates...")
    # this will call cnvs seperately if they end and start on the same coordinate, e.g. 1-1000 and 1000-4000 will be called as non-overlapping CNVs
    distinct_cnvs = pyranges.merge(slack=-1, by="type")
    return distinct_cnvs

def convert_to_df(cnv_pr):
    ''' convert the filtered pyranges object containing filtered
        CNVs back into a dataframe '''
    cnv_df = cnv_pr.as_df()
    # convert the type of chromsome column so we can merge
    cnv_df['Chromosome'] = cnv_df['Chromosome'].astype('str')
    # drop the var type column to make the pandas merge easier, add binSize column for traceability
    return cnv_df

def main():
    array_files = glob.glob("*_array_cnvs.csv")
    # file containing local_id\tfamily_id\tnhs_number. allows us to link the array report to the wgs sample
    sw_families = pd.read_csv('sw_families.tsv', sep="\t")

    for f in array_files:
        nhs = f.split('_')[0]
        grch38_bed = nhs+'_grch38.bed'
        # sometimes the table extraction from the pdf file doesn't work properly, this attempts to mitigate for that
        try:
            array_cnvs = pd.read_csv(f)
        except:
            array_cnvs = pd.read_csv(f, skiprows=3, names=['number', 'Chromosome', 'cytogenetic_location', 
                                'Start_grch37', 'Stop_grch37', 'Size', 'Type', 'Classification_initial', 'Classification_final'])
        # merge the grch37 array cnvs with the liftover grch38 cnvs based on the cnv id number. the output from liftover is a bed file missing key information
        grch38_cnvs = pd.read_csv(grch38_bed, sep="\t", header=None, names=['Chromosome', 'Start','End','number','chunks'])
        cnvs = array_cnvs.merge(grch38_cnvs, on="number", how='left', suffixes=('', '_y'))
        cnvs['Chromosome'] = cnvs['Chromosome'].astype('str')
        # find the ex number from nhs number
        ex_number = sw_families[sw_families['nhs'] == int(nhs)]['ex_number'].item()
        print(nhs + " " + ex_number)
        # read in fakealamut file
        fakealamut_file = glob.glob('*' + ex_number + '*savvycnv.fakealamut.txt')
        alamut_df = pd.read_csv(fakealamut_file[0], sep='\t')
        alamut_df['chrom'] = alamut_df['chrom'].astype('str')
        # convert array cnvs to pyranges object
        array_df = cnvs[['Chromosome','Start','End','number']]
        array_df = array_df.dropna()
        array_pr = pr.PyRanges(df=array_df)
        # convert fakealamut to pyranges object
        alamut_pr = alamut_to_pyranges(alamut_df)
        # join
        joined_pr = array_pr.join(alamut_pr, suffix="_savvycnv", report_overlap=True)
        if joined_pr:
            joined_df = convert_to_df(joined_pr)
            # merge with the alamut df to get savvycnv call information
            alamut_joined = alamut_df.merge(joined_df,
                                        left_on=['chrom', 'gDNAstart', 'gDNAend'],
                                        right_on=['Chromosome', 'Start_savvycnv', 'End_savvycnv'],
                                        how='inner', suffixes=('', '_y'))
            # remove duplicate rows (due to genes, transcripts etc)
            drop_columns = ['gene','geneId','transcript','exon','distNearestSS','alleleFreq','dgvEntry']
            alamut_joined_drop = alamut_joined.drop(columns=drop_columns).drop_duplicates()
            # join with the original grch38 array cnv file
            array_alamut_joined = cnvs.merge(alamut_joined_drop, left_on=['Chromosome', 'Start','End'], right_on=['chrom','Start','End'], how="left", suffixes=('', '_y'))
            array_alamut_joined.drop(array_alamut_joined.filter(regex='_y$').columns, axis=1, inplace=True)
            array_alamut_joined.rename({'Start': 'Start_grch38', 'End': 'Stop_grch38'}, axis=1, inplace=True)
            array_alamut_joined.to_csv("array_overlap_savvycnv/"+nhs+"_array_savvycnv.csv",index=False)
        else:
            cnvs.drop(cnvs.filter(regex='_y$').columns, axis=1, inplace=True)
            cnvs.rename({'Start': 'Start_grch38', 'End': 'Stop_grch38'}, axis=1, inplace=True)
            cnvs.to_csv(nhs+"_array_savvycnv.csv",index=False)
        # now get array cnvs that did not lift over
        cnvs.drop(cnvs.filter(regex='_y$').columns, axis=1, inplace=True)
        cnvs.rename({'Start': 'Start_grch38', 'End': 'Stop_grch38'}, axis=1, inplace=True)
        to_lift = cnvs[cnvs['Start_grch38'].isna()]
        if not to_lift.empty:
            to_lift.to_csv(nhs+"_array_to_liftover.csv",index=False)

if __name__=="__main__":
    main()