# Nupack-Analysis

This repository is meant to be a guide for using python to run multiple analyses on Nupack models.

![NUPACK_results_20230413_213625.png](figures%2FNUPACK_results_20230413_213625.png)
![NUPACK_results_20230413_213627.png](figures%2FNUPACK_results_20230413_213627.png)

## Project Overview
nupack-4.0.0.28

Nupack documentation: https://piercelab-caltech.github.io/nupack-docs/analysis/

If you want to know about the physical model and algorithm in this program, visit: http://www.nupack.org/home/model

This python code will run a Nupack analysis for many experimental conditions at once, and will produce a percent bound graph for the input conditions. The program will also generate dataframe of all analysis data that can be downloaded as a csv.

It is recommended to use the run_expt_01.py for general purpose use.

The code can accept a single 5 prime strand and can automatically generates a compliment strand.

To run the code, follow the instructions below for setting up the Docker container.

The run_expt_01.py script has the following input fields.

* Five prime strand sequence
* three prime strand sequence
* Nucleic acid type: DNA or RNA
* Temperature start in ˚C
* Temperature end in ˚C
* Temperature incriment
* Copy number of target/template in sample [Can be a list]
* Sample volume in μL [Can be a list]
* Sodium (Na+) molar concentration
* Magnesium (Mg++) molar concentration

Figure are automatically generated and saved in the figures directory.

Note: Users should input their own sequences and conditions, and run the new analysis.

## Getting Started
1. To get started, download or clone the repository. You can do this with by typing the following into terminal:

```bash
git clone https://github.com/dleon86/nupack_analysis.git
cd ~/nupack_analysis
```

2. Install Docker on your computer by following the instructions for your operating system on the Docker website: 
https://www.docker.com/get-started.

* Note: Docker for desktop needs to be running for the container to work.

3. Open a terminal or command prompt and pull the nptb image from Docker Hub by running the following command:

```bash
docker pull dleon86/nptb
```

5. Start a Docker container and mount the nupack_analysis directory to the container at /home/nupack_user/app. You can 
do this by running the following command:

```bash
docker run -it -v /path/to/nupack_analysis:/home/nupack_user/app -p 5000:5000 dleon86/nptb
```
Replace /path/to/nupack_analysis with the actual path to the nupack_analysis directory on your computer.

Once the container starts, you should see the terminal prompt change to `(base) 
root@xxxxxxxxxxxx:/home/nupack_user/app#` where the `x`'s are replaced with letters and numerals. 

* note: If you have previously created and started the container, you can reopen it with:

```bash
docker start -ai nptb
```

6. From here, you can run the experiment by running the following command, in the container:
```bash
python run_expt_01.py
```
This will execute the run_expt_01.py script, which will generate the figures in the figures directory.

If you want to run a Jupyter notebook, you can start the Jupyter server by running the following command in the container:

```css
jupyter notebook --ip=0.0.0.0 --port=5000 --no-browser
```
This will start the Jupyter server on port 5000 and allow you to access it from your local computer by going to http://localhost:5000 in your web browser.

7. To clean up the container when you are done using it, you can exit the container by running the exit command:
```bash
exit
```

or pressing Ctrl+D. This will stop the container and return you to your local terminal.

If you want to come back to the container later and continue working on the experiment, you can start the container 
again by running the same docker run command as before (step 5.b). This will mount the nupack_analysis directory to the 
container and allow you to continue working where you left off.

## Citation
If you use this code for generating data for publication, remember to cite the NUPACK authors with the citations below.

M.E. Fornace, N.J. Porubsky, and N.A. Pierce (2020). A unified dynamic programming framework for the analysis of interacting nucleic acid strands: enhanced models, scalability, and speed. ACS Synth Biol, 9:2665-2678, 2020. (pdf, supp info)

R. M. Dirks, J. S. Bois, J. M. Schaeffer, E. Winfree, and N. A. Pierce. Thermodynamic analysis of interacting nucleic acid strands. SIAM Rev, 49:65-88, 2007. (pdf)

## Contact
Daniel Leon - dleon@uw.edu GitHub: dleon86

or

Joe Henthorn - JosefH1@uw.edu GitHub: JoeHenthorn