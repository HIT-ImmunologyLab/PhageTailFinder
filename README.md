# PhageTailFinder: a software for phage tail module recognition based on hidden markov model
## Table of Contents
- [Introductiom](#introudction)
- [Requirements](#requirements)
- [Install](#install)
- [Usage](#usage)
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
```
git clone https://github.com/HIT-ImmunologyLab/PhageTailFinder
```

## Usage
### Command line options

![image](https://github.com/HIT-ImmunologyLab/PhageTailFinder/blob/main/image/useage.png)

### Start PhageTailFinder

The python script is also provided for expert users<br>
predict phage tail protein with default parameters:

```
predict.py [-h] -i INPUT -o OUTPUT [--accurate--mode] [--hmmscan--evalue HMMSCAN__EVALUE] [--phagelist--mode]

```
### Download HMM model to annotate tail protein or non-protein
wget http://www.microbiome-bigdata.com/PHISDetector/static/download/PhageTailFinder/db.tar.gz

### Outputs

File Name | Description
---|---
 
## Contributors
This project exists thanks to all the people who contributed.

## License

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
