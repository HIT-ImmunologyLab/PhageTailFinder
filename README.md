# PhageTailFinder: A software for phage tail module recognition based on hidden markov model
## Table of Contents
- [Introductiom](#introudction)
- [Requirements](#requirements)
- [Install](#install)
- [Usage](#usage)
- [Visualizations](#visualization)
- [Contributors](#contributors)
- [License](#license)
## Introudction
PhageTailFinders is an effective tool to predict phage tail protein based on hidden Markov model. The reliability of this tool has been tested, and the results have good accuracy and low error rate. At the same time, compared with other tools, the prediction speed of this tool is faster and the effect is better.The standalone version of the whole tool is built and good visualization results are provided. The Python program combining multiple steps is also concentrated in a script file, which is convenient for users to use the multi-phage batch prediction provided in the tool.
## Requirements ##
The source code is written by python3. <br>
First, please install the following python packages:

1. numpy
 
2. xml.etree.ElementTree
 
3. sklearn

4. matplotlib

## Install ##
#### Linux
- step1:Download the whole packages and partial profiles from [https://github.com/HIT-ImmunologyLab/PhageTailFinder](https://github.com/HIT-ImmunologyLab/PhageTailFinder)


## Usage
### Command line options

![image](https://github.com/HIT-ImmunologyLab/PhageTailFinder/blob/main/image/useage.png)

### Start DBSCAN-SWA

The python script is also provided for expert users<br>
1.predict prophages of query bacterium with default parameters:

```
python <path>/dbscan-swa.py --input <bac_path> --output <outdir> --prefix <prefix>
```
2. predict prophages of query bacterium and no phage annotation:
```
python <path>/dbscan-swa.py --input <bac_path> --output <outdir> --prefix <prefix> --add_annotation none
```
3. predict prophages of query bacterium and detect the bacterium-phage interaction between the query bacterium and query phage:
```
python <path>/dbscan-swa.py --input <bac_path> --output <outdir> --prefix <prefix> --add_annotation <phage_path>
```
### Outputs

File Name | Description
---|---
\<prefix\>\_DBSCAN-SWA\_prophage\_summary.txt | the tab-delimited table contains predicted prophage informations including prophage location, specific phage-related key words, CDS number, infecting virus species by a majority vote and att sites
\<prefix\>\_DBSCAN-SWA\_prophage.txt | this table not only contains the information in <prefix>\_DBSCAN-SWA\_prophage\_summary.txt but also contains the detailed information of prophage proteins and hit parameters between the prophage protein and hit uniprot virus protein
<prefix>\_DBSCAN-SWA\_prophage.fna| all predicted prophage Nucleotide sequences in FASTA format
<prefix>\_DBSCAN-SWA\_prophage.faa| all predicted prophage protein sequences in FASTA format
**Phage Annotation**| if add\_annotation!=none, the following files are in "prophage\_annotation" 
<prefix>\_prophage\_annotate\_phage.txt | the tab-delimited table contains the information of prophage-phage pairs with prophage\_homolog\_percent, prophage\_alignment\_identity and prophage\_alignment\_coverage
<prefix>\_prophage\_annotate\_phage.txt | the table contains the detailed information of bacterium-phage interactions including blastp and blastn results 
## Contributors
This project exists thanks to all the people who contribute.

## License

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
