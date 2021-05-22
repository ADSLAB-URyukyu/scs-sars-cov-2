# Supplementary Information for Article .....
This site provides the following supplementary information
 * Source Codes for Human SCS Analysis,
 * Source Codes for SARS-CoV2 SCS Analysis,
 * Supplementary Data

## Source Codes for Human SCS Analysis
You can use the Python code easily with your computer or the Google Colab & Google Drive.

### How to use the Human SCS Analysis Program with your computer
1. Download the source code and Protein Datasets
  * and then locate Human_SCS_Analysis.py in the current directory and Protein data as ./ncbi_dataset/protein.aa  
2. Run jupyter notebook at your current directory
3. Start to use the program by importing Human_SCS_Application  
````python:
    import  Human_SCS_Analysis as hscs  
    hscs.initializeFromScsDataset()  
    # You can set all data to use the application   
    hscs.menu()
    # You can see the command list
    # For example, to show the basic information of the dataset
    hscs.showBasicInformation()
````

### How to use Human SCS Analysis Program with Google Colab
1. Download Human_SCS_Application and Protein Datasets
  * locate Human_SCS_Application.py at the Google Drive directory and Protein datasets at ./ncbi_dataset/protein.aa 
2. Open a new notebook in the Google Colab and mount the Google Drive directory
3. Start to use the application by importing Human_SCS_Application  
````python:
    import  Human_SCS_Application as hscs  
    hscs.initializeFromScsDataset()  
    # You can set all data to use the application   
    hscs.menu()
    # You can see the command list
    # For example, to show the basic information of the dataset
    hscs.showBasicInformation()
````


## Source Codes for SARS-CoV-2 SCS Analysis    

Under construction

## Supplementary Data
We provide the following supplementary data as Excel files in the repository:
1. Self-nonself assignment and Nonself extraction,
2. Nonself clusters,
3. SARS-CoV-2 variant proteomes,
4. Spike mutations