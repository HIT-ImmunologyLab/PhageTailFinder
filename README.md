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

Second, please install the following tools:
## Usage
PhageTailFinder is a phage tail protein prediction tool. A standalone version has been developed for use in this project.The program supports single phage file input or multiple phage inputs.You can enter a phage sequence file in FastA or GenBank format, and the program will automatically convert the format and predict the return result.You can also enter the GenomeID list of phages, and the program will automatically identify and download all the phage sequences in the list and make batch predictions, and return a table to store the predictions of all phages.
### 
## Contributors
This project exists thanks to all the people who contribute.

## License

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
