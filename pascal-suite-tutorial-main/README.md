# PaScal Suite Tutorial
This is a hands-on tutorial that presents the basics of the Parallel Scalability Suite (PaScal Suite), a toolset for assessing the scalability trends of parallel programs.

## Description
PaScal Suite integrates two tools, the PaScal Analyzer and the PaScal Viewer. It simplifies the execution, measurement, and comparison of several executions of a parallel program, allowing the analysis of scalability trends in configuration environments with different amounts of processing elements and different workloads with visual elements that help the user to understand the program’s behavior and recognize scalability bottlenecks that may require in-depth optimization analysis. The toolset can be used to help the development of parallel programs that run on a single shared-memory computational node. Future versions will support distributed programs running on many nodes.

## Download the tutorial's files
Use the following command to clone this repo and download the tutorial's files.
```bash
git clone git@gitlab.com:lappsufrn/pascal-suite-tutorial.git
```

## Installing PaScal Analyzer
PaScal Analyzer need to be installed in the machine that the parallel program under analysis needs to run. The PaScal Releases repo provides a Linux x86_64 binary version of the Analyzer, which can be quickly installed locally with the following commands:
```bash
wget -c https://gitlab.com/lappsufrn/pascal-releases/-/archive/master/pascal-releases-master.zip
unzip pascal-releases-master.zip
cd pascal-releases-master/
source env.sh
```

To see a short help on how to use the PaScal Analyzer:
```bash
pascalanalyzer -h
```

## Usage
To compile the tutorial's programs from within the tutorial's directory:
```bash
cd pascal-suite-tutorial/
make
```

To submit a job to the NPAD's cluster using the nodes reserved to the tutorial:
```bash
sbatch pascalanalyzer_job.sh 
```

To view the status of your jobs:
```bash
squeue -u $USER
```

## Authors and acknowledgment
This tutorial was prepared by [Samuel Xavier-de-Souza](http://dca.ufrn.br/~samuel) with help from Anderson B. N. da Silva and Vitor R. G. da Silva.

## License
PaScal Suite Tutorial © 2023 by Samuel Xavier de Souza is licensed under [Attribution-NonCommercial-ShareAlike 4.0 International](http://creativecommons.org/licenses/by-nc-sa/4.0/).
