# PhageTailFinder: A software for phage tail module recognition based on hidden markov model
## Table of Contents
- [Background](#background)
- [Requirements](#requirements)
- [Install](#install)
- [Usage](#usage)
- [Visualizations](#visualization)
- [Contributors](#contributors)
- [License](#license)
## Background
PhageTailFinders is an effective tool to predict phage tail protein based on hidden Markov model.Hidden Markov model is a statistical probability model, which can be used to represent an observation sequence. The standalone version of the whole tool is built and good visualization results are provided. The Python program combining multiple steps is also concentrated in a script file, which is convenient for users to use the multi-phage batch prediction provided in the tool and can also provide more results and more accurate rules for researchers to use
## Requirements ##
The source code is written by python3. In addition, several tools have been applied in DBSCAN-SWA. Among these, Prokka requires installtion by users. <br>
First, please install the following python packages:

1. numpy
 
2. Biopython
 
3. sklearn

Second, please install the following tools:
1. Prokka in https://github.com/tseemann/prokka<br>
```
git clone https://github.com/tseemann/prokka.git
# install the dependencies:
sudo apt-get -y install bioperl libdatetime-perl libxml-simple-perl libdigest-md5-perl
# install perl package
sudo bash
export PERL_MM_USE_DEFAULT=1
export PERL_EXTUTILS_AUTOINSTALL="--defaultdeps"
perl -MCPAN -e 'install "XML::Simple"'
# install the prokka databases
prokka --setupdb
# test the installed prokka databases
prokka --listdb
```
**warning**: Prokka needs blast+ 2.8 or higher, so we provide the version of blast+ in bin directory, the users can install a latest blast+ and add it to the environment or use the blast+ provided by DBSCAN-SWA. Please ensure the usage of blast+ in your environment by eg: 
```
which makeblastdb
```



## License

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
