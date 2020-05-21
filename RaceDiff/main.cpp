//
//  main.cpp
//  RaceDiff
//
//  Created by Gang Peng on 5/21/20.
//  Copyright © 2020 Gang Peng. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <math.h>

#include "normal.h"

using namespace std;

int main(int argc, char ** argv) {
/*
 * -i input.txt // input file containing MAF information
 * -o output.txt // output file
 * -r race1 race2 race3 ... // race groups to compare. these groups should be in the first row (header) of the input file.
 * -c cutoff // significant cutoff. If the MAF difference is larger than this cutoff, the variant is significant differnt between the two race groups. default: 0.05
 */
    
    vector<string> mustOptions = {"-i", "-o", "-r"};
    vector<string> allOptions = {"-i", "-o", "-r", "-c"};
    
    double cutoff = 0.05;
    
    map<string, vector<string> > cmLine = parseCMLine(argc, argv, allOptions, mustOptions);
    
    string fnameIn = cmLine["-i"][0];
    string fnameOut = cmLine["-o"][0];
    
    ifstream fin(fnameIn.c_str());
    if(!fin.is_open()){
        cout << "Cannot open file " << fnameIn << ". Please check the input file name." << endl;
        return -1;
    }
    
    ofstream fout(fnameOut.c_str());
    if(!fout.is_open()){
        cout << "Cannot open file " << fnameOut << ". Please check the output file name." << endl;
        return -1;
    }
    
    string header;
    getline(fin, header);
    
    //check header and race in the header
    vector<string> vsHeader = split(header, "\t");
    
    vector<size_t> idxRace;
    for(size_t i=0; i<cmLine["-r"].size(); i++){
        vector<string>::iterator it = find(vsHeader.begin(), vsHeader.end(), cmLine["-r"][i]);
        if(it!=vsHeader.end()){
            idxRace.push_back(distance(vsHeader.begin(), it));
        } else {
            cout << "Cannot find race \"" << cmLine["-r"][i] << "\" in the header of the input file. Plese take a look at the header of the input file." << endl;
            return -1;
        }
    }
    
    size_t idxSymbol=0;
    vector<string>::iterator itSymbol = find(vsHeader.begin(), vsHeader.end(), "SYMBOL");
    if(itSymbol == vsHeader.end()){
        cout << "Cannot find SYMBOL column in input file. Please make sure \"" << fnameIn <<"\" is the right input file." <<endl;
        return -1;
    } else {
        idxSymbol = distance(vsHeader.begin(), itSymbol);
    }
    
    if(idxRace.size() == 1){
        fout<<"Chr\tPos\tSymbol";
        fout<<"\t"<<cmLine["-r"][0]<<"-Other"<<endl;
        
        string currentGene = "";
        string currentChr = "";
        vector<int> pos;
        
        int numSigPositive = 0;
        int numSigNegative = 0;
        int numTotal = 0;
        
        int numSigPositiveTotal = 0;
        int numSigNegativeTotal = 0;
        int numTotalTotal = 0;
        
        
    } else {
        //fout header of output file
        size_t numCompare = (cmLine["-r"].size()-1) * (cmLine["-r"].size()) / 2;
        fout<<"Chr\tPos\tSymbol";
        for(size_t i=0; i < (cmLine["-r"].size()-1); i++){
            for(size_t j=(i+1); j<cmLine["-r"].size(); j++){
                fout<<"\t"<<cmLine["-r"][i]<<"-"<<cmLine["-r"][j];
            }
        }
        fout<<endl;
        
        string currentGene = "";
        string currentChr = "";
        vector<int> pos;
        
        vector<int> numSigPositive(numCompare, 0);
        vector<int> numSigNegative(numCompare, 0);
        vector<int> numTotal(numCompare, 0);
        
        vector<int> numSigPositiveTotal(numCompare, 0);
        vector<int> numSigNegativeTotal(numCompare, 0);
        vector<int> numTotalTotal(numCompare, 0);
        
        while(!fin.eof()){
            string fline;
            getline(fin, fline);
            
            if(fline.size() < 2){
                break;
            }
            
            vector<string> vsLine = split(fline, "\t");
            if(vsLine[idxSymbol]==""){
                continue;
            }
            
            if(vsLine[idxSymbol] != currentGene){
                if(currentGene!=""){
                    fout << currentChr;
                    if(pos.size() % 2 == 0){
                        int posTmp = (int)(pos[pos.size()/2] + pos[pos.size()/2-1])/2;
                        fout << "\t" << posTmp << "\t" << currentGene;
                    } else {
                        int posTmp = pos[pos.size()/2];
                        fout << "\t" << posTmp << "\t" << currentGene;
                    }
                    
                    for(size_t i=0; i<numCompare; i++){
                        fout << "\t" << numSigPositive[i] << "|" << numSigNegative[i] << "|" << numTotal[i];
                    }
                    fout <<endl;
                }
                currentGene = vsLine[idxSymbol];
                currentChr = vsLine[0];
                pos.clear();
                fill(numSigPositive.begin(), numSigPositive.end(), 0);
                fill(numSigNegative.begin(), numSigNegative.end(), 0);
                fill(numTotal.begin(), numTotal.end(), 0);
            }
            
            
            pos.push_back(stoi(vsLine[1]));
            int idx = 0;
            for(size_t i=0; i<(idxRace.size()-1); i++){
                for(size_t j=(i+1); j<idxRace.size(); j++){
                    if(vsLine[idxRace[i]] != "NA" && vsLine[idxRace[i]] != "NA"){
                        double maf1 = stof(vsLine[idxRace[i]]);
                        double maf2 = stof(vsLine[idxRace[j]]);
                        if(fabs(maf1-maf2) > cutoff){
                            if(maf1 > maf2){
                                numSigPositive[idx]++;
                                numSigPositiveTotal[idx]++;
                            } else {
                                numSigNegative[idx]++;
                                numSigNegativeTotal[idx]++;
                            }
                        }
                        numTotal[idx]++;
                        numTotalTotal[idx]++;
                    }
                    idx++;
                }
            }
        }
        
        //output last record
        fout << currentChr;
        if(pos.size() % 2 == 0){
            int posTmp = (int)(pos[pos.size()/2] + pos[pos.size()/2-1])/2;
            fout << "\t" << posTmp << "\t" << currentGene;
        } else {
            int posTmp = pos[pos.size()/2];
            fout << "\t" << posTmp << "\t" << currentGene;
        }
        
        for(size_t i=0; i<numCompare; i++){
            fout << "\t" << numSigPositive[i] << "|" << numSigNegative[i] << "|" << numTotal[i];
        }
        fout <<endl;
        
        //output total
        fout<<"Total\t0\tTotal";
        for(size_t i=0; i<numCompare; i++){
            fout << "\t" << numSigPositiveTotal[i] << "|" << numSigNegativeTotal[i] << "|" << numTotalTotal[i];
        }
        fout <<endl;
    }
    
    fin.close();
    fout.close();
    
    return 0;
}
