# Haplotype Phasing

A collection of algorithms that take genotypes of individuals coming from admixed populations as input and outputs their haplotypes.

## Usage
### pip3
If you don't already have pip3 installed, you'll need to run the following command in your terminal before proceeding:
```
$ sudo apt-get install python3-pip
```
### Python dependencies
You'll need to install the following python modules:
```
$ pip3 install numpy
$ pip3 install progress
```
## Testing
### Clark's Algorithm
When you're inside the src folder, type
```
$ python3 clarks.py
```
to test on example data 1, with default block size = 180 and window length = 30.

### Switch Accuracy
Running the command
```
$ Rscript calculate_switch_accuracy.R ../data/example_data_1/example_data_1_my_sol.txt ../data/example_data_1/example_data_1_sol.txt
```
should output:
```
[1] "Calculating switch accuracy"
[1] "Switch accuray is:"
[1] 0.8472889
[1] "Thank you, goodbye"
```
