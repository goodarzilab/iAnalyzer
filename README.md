# iAnalyzeR
An analytical pipeline for powered analysis of CRISPR screens

### Setting up the environment
Use the following command to create a conda environment with the essential packages used by iAnalyzeR.
```bash
conda env create --file=iAnalyzeR.yaml
conda activate iAnalyzeR
```
### Run iAnalyzeR
#### Create the metadata file
The metadata file resides in the working directory and lists the required information for each sample (e.g. sample name, path to files, sample type, and *etc*). For example:

|---|---|---|| Sample.name  |    fastq         |lib.type|sample.type| sample.rep ||---|---|---|| hi_r1        | hi_r1.fastq.gz   |   F    |     high   |      1 || lo_r1        | ro_l1.fastq.gz   |   F    |     low    |      1 |

#### Run the analysis
The following command will then run the analysis:
```bash
python iAnalyzer.py --runMode --ref=<ref> <metadata> <formula> <outfile.txt>
```
#### Options
Run `python iAnalyzer.py` for usage.
The following are the options:
1. `--runMode` vs. `--printMode`: `--printMode` prints all the commands that are run
2. `-a` or `--aligner`: The aligner used (currently only bowtie2 supported)
3. `-w` or `--weighting`: Method for combining z-score. Options are 'n' (sample size), 'SE' (standard error), and 'SES' (standardized effect size) (default SES)
4. `--hasUMI` or `--noUMI`: Whether reads contain UMI or not
5. `--reflib`: The bowtie2 index basename (e.g. CRISPRi_v2_human_library for CRISPRi_v2_human_library.fa)
6. `--ref`: The reference samples that others are compare to (e.g. 'low' in example above)
7. `formula`: The design formula (e.g. ~sample.type in example above, i.e. high vs. low)
