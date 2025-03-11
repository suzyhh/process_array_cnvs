# Comparing CNVs detected by array to CNVs detected by SavvyCNV

## Set up
Developed using python 3.9
```
python -m venv my_venv
source my_venv/bin/activate
pip install -r requirements.txt
```

Also required is a local installation of UCSC's [liftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver) tool along with your desired chain file

## Running
There are three scripts involved:
* `get_pdf_tables.py`: Reads PDF array report files from NBT and extracts the table containing detected CNVs.
* `merge_builds.py`: Takes lifted over array CNVs and intersects these with CNVs detected for that sample by WGS using SavvyCNV to find overlapping CNVs
* `liftOver_cnvs.sh`: Control script that runs the two previously mentioned python scripts and also performs liftover of the array CNVs which are in GRCh37 to GRCh38.

You will also need to have your SavvyCNV files available. The script is set up to take SavvyCNV files annotated with Alamut Batch
