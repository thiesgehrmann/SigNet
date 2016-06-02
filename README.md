# SigNet
A unified testing and selection method to discover statistically significantly mutated sub-networks

### Usage

The pipeline.sh file contains the four parts of this tool that need to be executed in that order.
 1. download_data.sh # Download the data
 2. preprocessing.py # Perform preprocessing of the data
 3. diffusion.m      # Perform diffusion over the network
 4. analysis.py      # Discover significant clusters and downstream analysis

### Dependencies

* gffread utility (http://cole-trapnell-lab.github.io/cufflinks/file_formats/#the-gffread-utility)
* fastcluster (https://pypi.python.org/pypi/fastcluster)
* Ibidas (https://pypi.python.org/pypi/Ibidas)
* Numpy
* Scipy
* DTcut (https://github.com/thiesgehrmann/DTcut)
* Python 2.7
* matplotlib
* seaborn
* pandas
* bx-python (https://pypi.python.org/pypi/bx-python)


