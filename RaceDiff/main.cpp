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
 * -r race1 race2 race3 ... // race groups to compare. these groups should be in the first row (header) of the input file. -r all: comparing all pairs  -r other: comparing every race to other
 * -c cutoff // significant cutoff. If the MAF difference is larger than this cutoff, the variant is significant differnt between the two race groups. default: 0.05
 */
    
    vector<string> mustOptions = {"-i", "-o", "-r"};
    vector<string> allOptions = {"-i", "-o", "-r", "-c"};
    
    double cutoff = 0.05;
    
    map<string, vector<string> > cmLine = parseCMLine(argc, argv, allOptions, mustOptions);
    if(cmLine.size()==0){
        return -1;
    }
    
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
    if(cmLine["-r"].size()==1){
        if(cmLine["-r"][0] == "other"){
            idxRace.clear();
        } else if(cmLine["-r"][0] == "all"){
            idxRace.clear();
            // There is a bug here
            // I assume race info from 7 to 11 // 6-10
            for(size_t i=6; i<11; i++){
                idxRace.push_back(i);
            }
        } else {
            vector<string>::iterator it = find(vsHeader.begin(), vsHeader.end(), cmLine["-r"][0]);
            if(it!=vsHeader.end()){
                idxRace.push_back(distance(vsHeader.begin(), it));
            } else {
                cout << "Cannot find race \"" << cmLine["-r"][0] << "\" in the header of the input file. Plese take a look at the header of the input file." << endl;
                return -1;
            }
        }
    } else {
        for(size_t i=0; i<cmLine["-r"].size(); i++){
            vector<string>::iterator it = find(vsHeader.begin(), vsHeader.end(), cmLine["-r"][i]);
            if(it!=vsHeader.end()){
                idxRace.push_back(distance(vsHeader.begin(), it));
            } else {
                cout << "Cannot find race \"" << cmLine["-r"][i] << "\" in the header of the input file. Plese take a look at the header of the input file." << endl;
                return -1;
            }
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
    
    cout << "Start" <<endl;
    
    if(idxRace.size() == 0){
        // There is a bug here
        // I assume race info from 7 to 11 // 6-10
        int numCompare = 5;
        fout<<"Chr\tPos\tSymbol";
        for(size_t i=6; i<11; i++){
            fout << "\t" << vsHeader[i] << "-Other";
        }
        fout << endl;
        
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
            for(size_t i=6; i<11; i++){
                if(vsLine[i] == "NA"){
                    continue;
                }
                double maf1 = stof(vsLine[i]);
                double maf2 = 0;
                int numOther = 0;
                for(size_t j=6; j<11; j++){
                    if(j!=i){
                        if(vsLine[j] != "NA"){
                            maf2 += stof(vsLine[j]);
                            numOther++;
                        }
                    }
                }
                if(numOther==0){
                    continue;
                }
                
                maf2 = maf2/numOther;
                if(fabs(maf1-maf2) > cutoff){
                    if(maf1 > maf2){
                        numSigPositive[i-6]++;
                        numSigPositiveTotal[i-6]++;
                    } else {
                        numSigNegative[i-6]++;
                        numSigNegativeTotal[i-6]++;
                    }
                }
                numTotal[i-6]++;
                numTotalTotal[i-6]++;
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
        
    } else if(idxRace.size() == 1){
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
                    
                    fout << "\t" << numSigPositive << "|" << numSigNegative << "|" << numTotal << endl;
                }
                currentGene = vsLine[idxSymbol];
                currentChr = vsLine[0];
                pos.clear();
                numSigPositive = 0;
                numSigNegative = 0;
                numTotal = 0;
            }
            
            
            pos.push_back(stoi(vsLine[1]));
            // There is a bug here
            // I assume race info from 7 to 11 // 6-10
            if(vsLine[idxRace[0]] == "NA"){
                continue;
            }
            
            double maf1 = stof(vsLine[idxRace[0]]);
            
            double maf2 = 0;
            int numOther = 0;
            for(size_t i=6; i<11; i++){
                if(i != idxRace[0]){
                    if(vsLine[i] != "NA"){
                        maf2 += stof(vsLine[i]);
                        numOther++;
                    }
                }
            }
            
            if(numOther==0){
                continue;
            }
            
            maf2 = maf2/numOther;
            
            if(fabs(maf1-maf2) > cutoff){
                if(maf1 > maf2){
                    numSigPositive++;
                    numSigPositiveTotal++;
                } else {
                    numSigNegative++;
                    numSigNegativeTotal++;
                }
            }
            numTotal++;
            numTotalTotal++;
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
        
        fout << "\t" << numSigPositive << "|" << numSigNegative << "|" << numTotal << endl;
        
        //output total
        fout << "Total\t0\tTotal\t" << numSigPositiveTotal << "|" << numSigNegativeTotal << "|" << numTotalTotal << endl;
        
    } else {
        //fout header of output file
        size_t numCompare = (idxRace.size()-1) * (idxRace.size()) / 2;
        fout<<"Chr\tPos\tSymbol";
        for(size_t i=0; i<(idxRace.size()-1); i++){
            for(size_t j=(i+1); j<idxRace.size(); j++){
                fout<<"\t"<<vsHeader[idxRace[i]]<<"-"<<vsHeader[idxRace[j]];
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
                    if(vsLine[idxRace[i]] != "NA" && vsLine[idxRace[j]] != "NA"){
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
    
    cout << "Finished!" << endl;
    
    return 0;
}
