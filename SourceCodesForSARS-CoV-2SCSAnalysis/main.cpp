#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>

#define HUMAN_FILE "./ncbi_dataset/protein.faa" //Path of the human file 
#define COVID_FILE "./ncbi_dataset/coronavirus" //Path of the covid-19 file
#define SCS3 10648   
#define SCS4 234256  
#define SCS5 5153632 
#define SIZE_AVALIABILITIES3(amino3_avaliabilities) (sizeof(scs3_avaliabilities)/sizeof(scs3_avaliabilities[0]))
#define SIZE_AVALIABILITIES4(amino4_avaliabilities) (sizeof(scs4_avaliabilities)/sizeof(scs4_avaliabilities[0]))
#define SIZE_AVALIABILITIES5(amino5_avaliabilities) (sizeof(scs5_avaliabilities)/sizeof(scs5_avaliabilities[0]))
#define SIZE_OCCURRENCES3(amino3_occurrences) (sizeof(scs3_occurrences)/sizeof(scs3_occurrences[0]))
#define SIZE_OCCURRENCES4(amino4_occurrences) (sizeof(scs4_occurrences)/sizeof(scs3_occurrences[0]))
#define SIZE_OCCURRENCES5(amino5_occurrences) (sizeof(scs5_occurrences)/sizeof(scs5_occurrences[0]))
using namespace std;

vector<string> split_naive(const string &s, char delim);

extern void calculation_occurrence(int length, int* occurrences, string file);                                                       
extern void extraction_zero_occurrences(int length, int* scs_occurrences, string file);                                                                                         

//g++ -std=c++17 -stdlib=libc++ -o main main.cpp calculation.cpp transform.cpp 
int main(){
    static int amino_occurrences[22];         
    static int scs3_occurrences[SCS3];      
    static int scs4_occurrences[SCS4];       
    static int scs5_occurrences[SCS5];       
    static int element_occurrences3[SCS3];     
    static int element_occurrences4[SCS4];     
    static int element_occurrences5[SCS5];     
    static double scs3_avaliabilities[SCS3]; 
    static double scs4_avaliabilities[SCS4]; 
    static double scs5_avaliabilities[SCS5]; 
    static int element_avaliabilities3[SCS3];  
    static int element_avaliabilities4[SCS4];  
    static int element_avaliabilities5[SCS5];  
    std::ifstream ifs(HUMAN_FILE);              
    
    if (ifs.fail()) { 
        std::cerr << "Failed to open file." << std::endl;
        return -1;
    }
    for (int i = 0; i < SCS5; i++){ 
        if (i < SCS4){
            if (i < SCS3){
                if (i < 22){
                    amino_occurrences[i] = 0;
                }
                scs3_occurrences[i] = 0;
                element_occurrences3[i] = i;
                scs3_avaliabilities[i] = 0;
                element_avaliabilities3[i] = i;
            }
            scs4_occurrences[i] = 0;
            element_occurrences4[i] = i;
            scs4_avaliabilities[i] = 0;
            element_avaliabilities4[i] = i;
            }
        scs5_occurrences[i] = 0;
        element_occurrences5[i] = i;
        scs5_avaliabilities[i] = 0;
        element_avaliabilities5[i] = i;
    }
   
    calculation_occurrence(5, scs5_occurrences, HUMAN_FILE); 
    for (const auto & file : std::filesystem::directory_iterator(COVID_FILE)){
        extraction_zero_occurrences(5, scs5_occurrences, file.path());    
    }
    return 0;
}
