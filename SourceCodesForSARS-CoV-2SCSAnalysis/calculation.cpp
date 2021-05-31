#include <stdio.h>
#include <iostream>
#include <cmath>
#include <string>
#include <tuple>
#include <fstream>
#include <sstream>
#include <vector> 
using namespace std;

extern int check_amino(char amino);
extern int transform_scs_10(int n, string amino5);         
extern std::string transform_10_scs(int element, int n);

vector<string> split_naive(const string &s, char delim) {
    vector<string> elems;
    string item;
    for (char ch: s) {
        if (ch == delim) {
            if (!item.empty())
                elems.push_back(item);
            item.clear();
        }
        else {
            item += ch;
        }
    }
    if (!item.empty())
        elems.push_back(item);
    return elems;
}

void calculation_occurrence(int length, int* occurrences, string file){             
    int protein_num = 0;                  
    std::string scs;                   
    std::string protein;                   
    std::string str;                       
    std::ifstream ifs(file);
    int ddd = 0;             
    
    while (getline(ifs, str)) { 
        if (str[0] == '>'){     
            protein_num ++;

            if (protein_num != 1){
                for (int i = 0; i < protein.length()-length+1; i++){
                    scs = protein[i];
                    for (int i2 = 1; i2 < length; i2++){                        
                        scs += protein[i+i2];
                    }
                    for (int j = 0; j < 5; j++){    
                        if (scs[j] == 'X' || scs[j] == 'U'){
                            ddd ++;
                            break;
                        }              
                    }
                    occurrences[transform_scs_10(length, scs)] ++;
                }
            }
            protein = "";
        }else{
            protein += str;
        }
    }
    for (int i = 0; i < protein.length()-length+1; i++){
        scs = protein[i];
        for (int i2 = 1; i2 < length; i2++){                        
            scs += protein[i+i2];
        }
        occurrences[transform_scs_10(length, scs)] ++;
    }
}
void extraction_zero_occurrences(int length, int* scs_occurrences, string file){
    int protein_num = 0;       
    int index = 0;       
    static int zero_covid19[5153632];
    std::string protein;   
    std::string scs;          
    std::string str;                       
    std::ifstream ifs(file);    
    vector<std::string> strvec = split_naive(file, '_');
    ofstream outputfile("./SARS-CoV-2_SCS_Analysis/"+strvec[3]+".csv");
    for (int i = 0; i < std::pow(22, length); i++){
        zero_covid19[i] = 0;
    }
    while (getline(ifs, str)) { 
        if (str[0] == '>'){  
            protein_num ++;
            if (protein_num != 1){
                for (int i = 0; i < protein.length() - length + 1; i++){
                    scs = protein[i];
                    
                    for (int i2 = 1; i2 < length; i2++){                        
                        scs += protein[i+i2];
                    }
                    if (scs_occurrences[transform_scs_10(length, scs)] == 0){
                        zero_covid19[transform_scs_10(length, scs)] ++;
                        outputfile<<protein_num-1;
                        outputfile<<"-";
                        outputfile<<i+1;
                        outputfile<<",";
                        outputfile<<1;
                        outputfile<<",";
                        outputfile<<scs[0]<<endl;
                    }else{
                        outputfile<<protein_num-1;
                        outputfile<<"-";
                        outputfile<<i+1;
                        outputfile<<",";
                        outputfile<<0;
                        outputfile<<",";
                        outputfile<<scs[0]<<endl;
                    }
                }
                protein = "";
            }
        }else{
            protein += str;
        }
    }
    for (int i = 0; i < protein.length() - length + 1; i++){
        scs = protein[i];
        
        for (int i2 = 1; i2 < length; i2++){                        
            scs += protein[i+i2];
        }
        if (scs_occurrences[transform_scs_10(length, scs)] == 0){
            zero_covid19[transform_scs_10(length, scs)] ++;
            outputfile<<protein_num;
            outputfile<<"-";
            outputfile<<i+1;
            outputfile<<",";
            outputfile<<1;
            outputfile<<",";
            outputfile<<scs[0]<<endl;
        }else{
            outputfile<<protein_num;
            outputfile<<"-";
            outputfile<<i+1;
            outputfile<<",";
            outputfile<<0;
            outputfile<<",";
            outputfile<<scs[0]<<endl;;
        }
    }
}

