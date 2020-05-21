/*
 * normal.cpp
 *
 *  Created on: Jun 21, 2011
 *      Author: gpeng
 */

#include "normal.h"

#include <iostream>

using namespace std;

std::vector<std::string> split(const std::string & src, const std::string & tok, bool trim, std::string NulStr)
{
	std::vector <std::string> rlt;
	if(src.empty() || tok.empty())
	{
		return rlt;
	}


	size_t head=0, tail=0,length=0;
	while((tail=(src.find(tok,head))) != std::string::npos)
	{
		if((length=tail-head)>0)
		{
			rlt.push_back(src.substr(head,length));
		}
		else if(trim==false)
		{
			rlt.push_back(NulStr);
		}

		head=tail+1;
	}

	std::string endStr=src.substr(head);
	if(endStr.empty())
	{
		if(trim==false)
		{
			rlt.push_back(NulStr);
		}
	}
	else
	{
		rlt.push_back(endStr);
	}

	return rlt;
}


std::vector<std::string> split2(const std::string & src, const std::string & toks, bool trim, std::string NulStr)
{

	std::vector<std::string> rlt;
	if(src.empty() || toks.empty())
	{
		return  rlt;
	}


	size_t head=0, tail=0,length=0;
	while((tail=(src.find_first_of(toks,head))) != std::string::npos)
	{
		if((length=tail-head)>0)
		{
			rlt.push_back(src.substr(head,length));
		}
		else if(trim==false)
		{
			rlt.push_back(NulStr);
		}

		head=tail+1;
	}

	std::string endStr=src.substr(head);
	if(endStr.empty())
	{
		if(trim==false)
		{
			rlt.push_back(NulStr);
		}
	}
	else
	{
		rlt.push_back(endStr);
	}

	return rlt;
}

std::map<std::string, std::vector<std::string> > parseCMLine(int argc, char* argv[]) {
    map<string, vector<string> > rlt;
    for(int i=1; i<argc; i++){
        string para = argv[i];
        if(para[0] == '-'){
            vector<string> tmp;
            rlt[para] = tmp;
            while(true){
                i++;
                if(i >= argc){
                    break;
                }
                if(argv[i][0] == '-'){
                    i--;
                    break;
                }
                rlt[para].push_back(string(argv[i]));
            }
            if(i >= argc){
                break;
            }
        } else {
            cout << para << " is a parameter withou option and is discarded" << endl;
        }
    }
    return rlt;
}


std::map<std::string, std::vector<std::string> > parseCMLine(int argc, char* argv[], const std::vector<std::string> & allOptions, const std::vector<std::string> & mustOptions){
    map<string, vector<string> > rlt;
    vector<bool> tag;
    if(mustOptions.size() > 0){
        tag = vector<bool> (mustOptions.size(), false);
    }
    
    
    for(int i=1; i<argc; i++){
        string para = argv[i];
        if(para[0] == '-'){
            
            if(find(allOptions.begin(), allOptions.end(), para) != allOptions.end()){
                vector<string> tmp;
                rlt[para] = tmp;
                
                if(mustOptions.size() > 0){
                    for(size_t i=0; i<mustOptions.size(); i++){
                        if(mustOptions[i] == para){
                            tag[i] = true;
                        }
                    }
                }
                
                while(true){
                    i++;
                    if(i >= argc){
                        break;
                    }
                    if(argv[i][0] == '-'){
                        i--;
                        break;
                    }
                    rlt[para].push_back(string(argv[i]));
                }
            } else {
                cout << "Option " << para << " is not included in this program. This option and the following parameters after this option are discarded: " << endl;
                bool idxFirst = true;
                while(true){
                    i++;
                    if(i >= argc){
                        break;
                    }
                    if(argv[i][0] == '-'){
                        i--;
                        break;
                    }
                    if(idxFirst){
                        cout << argv[i];
                        idxFirst = false;
                    } else {
                        cout << ", " << argv[i];
                    }
                    
                }
                cout << endl;
            }
            
            if(i >= argc){
                break;
            }
        } else {
            cout << para << " is a parameter withou option and is discarded" << endl;
        }
    }
    
    bool mustOptionsNotSet = false;
    if(mustOptions.size() > 0){
        for(size_t i=0; i<tag.size(); i++){
            if(!tag[i]){
                cout << mustOptions[i] << " must be set in this program. " << endl;
                mustOptionsNotSet = true;
            }
        }
    }
    
    if(mustOptionsNotSet){
        rlt.clear();
    }
    
    return rlt;
}
