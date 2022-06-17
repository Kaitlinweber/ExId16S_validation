<div align="center">
    <h1>ExId16S validation</h1>
    <br />
    <h2>Script to validate the python package ExId16S</h2>
    <br />
</div>

## General information
* **Author:** Kaitlin Weber
* **Commissioned by:** Rijksinstituut voor Volksgezondheid en Milieu (RIVM)

## About this script

This python script validates the results of [ExId16S](https://github.com/Kaitlinweber/exid16s), by comparing the summary file kreport with the results of the gold standards 16S rDNA sequencing and WGS sequencing. The script creates bar plots, for comparison with the gold standards. In addition, the script also has the option to compare the results of the [custom made database](https://github.com/Kaitlinweber/combine_16S_database) with a standard SILVA database.

## Installation

1. Clone the repository:

```
git clone https://github.com/Kaitlinweber/exid16s_validation
```

2. Enter the directory with the pipeline and install the conda environment:

```
cd exid16s_validation
conda env install -f envs/vaidation_env.yaml
```

### Required parameters

* ```-i, --input```  Kraken summary kreport from ExId16S results, based on the custom made database.
* ```-w, --wgs``` Kraken summary kreport obtained from Kraken2 and Bracken with WGS data.
* ```-s, --sample_list``` Excel sheet with 16S rDNA Sanger sequencing results, with the sample number, genus name and species name in a seperate column. 
* ```-o, --output``` Pathway to output directory, if directory does not exists, directory will be created


### Optional parameters

* ```-si, --silva``` For the comparison of the created database with the SILVA database, the pathway to the SILVA database Kraken kreports can be given here


### The base command to run this script 

```
python validation.py -i [path/to/input/dir/ExId16S] -w [path/to/input/dir/WGS] -s [path/to/input/dir/16S_sequencing] -o [path/to/output/dir] 
```
