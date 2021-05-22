# Supplementary Information for Article .....
This site provides the following supplementary information
 * Human SCS Application,
 * SARS-CoV2 SCS Application,
 * Supplementary Data

## Human SCS Application
You can use the application easily with your computer or the Google Colaboratory & the Google Drive.

### How to use Human SCS Application with your computer
1. Download Human_SCS_Application and Protein Datasets
  * and then locate Human_SCS_Application.py in the current directory and Protein data as ./ncbi_dataset/protein.aa  
2. Run jupyter notebook at your current directory
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

### How to use Human SCS Application with Google Colaboratory
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


## SARS-CoV-2 SCS Application    

Under construction

## Supplementary Data
We provide the following supplementary data as Excel files in the repository:
1. Self-nonself assignment and Nonself extraction,
2. Nonself clusters,
3. SARS-CoV-2 variant proteomes,
4. Spike mutations
