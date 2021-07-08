
## About
This site is for Supplementary Information for Article *Self and nonself short constituent sequences of amino acids in the SARS-CoV-2 proteome for vaccine development* by Joji M. Otaki<sup>1</sup>, Wataru Nakasone<sup>2</sup>, Morikazu Nakamura<sup>2</sup>

1 The BCPH Unit of Molecular Physiology, Department of Chemistry, Biology and Marine Science, University of the Ryukyus, Okinawa 903-0213, Japan

2 Computer Science and Intelligent Systems Unit, Department of Engineering, Faculty of Engineering, University of the Ryukyus, Okinawa 903-0213, Japan


This site provides the following supplementary information
 * Source Codes for Human SCS Analysis,
 * Source Codes for SARS-CoV2 SCS Analysis,
 * Supplementary Data


## Source Codes for Human SCS Analysis
You can use the Python code easily with your computer or the Google Colab & Google Drive.

### How to use the Human SCS Analysis Program with your computer
1. Download the source code and Protein Datasets
  * and then locate Human_SCS_Analysis.py in the current directory and Protein data as ./ncbi_dataset/protein.faa  
2. Run jupyter notebook at your current directory
3. Start to use the program by importing Human_SCS_Application  
````python:
    import  Human_SCS_Analysis as hscs  
    hscs.initializeFromProteinDataset() 
    # You can set all data to use the application   
    hscs.menu()
    # You can see the command list
    # For example, to show the basic information of the dataset
    hscs.showBasicInformation()
````

### How to use Human SCS Analysis Program with Google Colab
1. Download the source code and Protein Datasets
  * locate Human_SCS_Analysis.py at the Google Drive directory and Protein datasets at ./ncbi_dataset/protein.faa 
2. Open a new notebook in the Google Colab and mount the Google Drive directory
3. Start to use the application by importing Human_SCS_Application  
````python:
    import  Human_SCS_Analysis as hscs   
    hscs.initializeFromProteinDataset()
    # You can set all data to use the application   
    hscs.menu()
    # You can see the command list
    # For example, to show the basic information of the dataset
    hscs.showBasicInformation()
````


## Source Codes for SARS-CoV-2 SCS Analysis 
0. Setup an environment on your computer to compile c++ source codes. Commandline tools are required. For linux and MacOS, the standard terminal tool is enough. For Windows 10 users, we recommend installing WSL or WSL2. 
1. Download SARS-CoV-2_SCS_Application and Protein Datasets to a working directory.
2. Run _make_ command at the working directory to build the software:
````console
    yourPC:~$ make
````
3. Instead of the _make_ command, you can compile the source codes directly to bilud the software:
````console
    yourPC:~$ g++ -std=c++17 -stdlib=libc++ -o covid-scs-analysis main.cpp calculation.cpp transform.cpp
````
4. You can run the application software with input data located at the directory ./ncbi_dataset by the following command:
````console
    yourPC:~$ ./covid-scs-analysis
````
5. The output files are generated under the output directory ./SARS-CoV-2_SCS_Analysis as csv files.

## Supplementary Data
We provide the following supplementary data as Excel files: [Supplementary Data](https://github.com/ADSLAB-URyukyu/scs-sars-cov-2/tree/main/Supplementary%20Data):
1. Self-nonself assignment and Nonself extraction,
2. Nonself clusters,
3. SARS-CoV-2 variant proteomes,
4. Spike mutations