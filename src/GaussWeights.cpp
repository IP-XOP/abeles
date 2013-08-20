//
//  GaussWeights.cpp
//  
//
//  Created by Andrew Nelson on 30/01/13.
//
//

#include "GaussWeights.h"
#include "gauss_legendre.h"

GaussWeight::GaussWeight(){};

GaussWeight::~GaussWeight(){};

GaussWeights::GaussWeights(){};

GaussWeights::~GaussWeights(){};

void GaussWeights::getGaussWeight(int n, std::vector<double> &x, std::vector<double> &w){
    std::map<int, GaussWeight>::iterator it;
    GaussWeight gw;
    int i, m;
    extern pthread_mutex_t changeWeightsMutex;

    pthread_mutex_lock(&changeWeightsMutex);
    it = gaussWeightsStorage.find(n);
    
    if(it == gaussWeightsStorage.end()){
        //get new gauss weights
        m = (n + 1) >> 1;
        gw.w.resize(m);
        gw.x.resize(m);
       
        gauss_legendre_tbl(n, &gw.x[0], &gw.w[0], 1.0e-17);
        
        std::reverse(gw.x.begin(),gw.x.end());
        std::reverse(gw.w.begin(),gw.w.end());
        
        for(i = m - 1 - (n & 1); i >= 0; i--){
            gw.x.push_back(-1. * gw.x[i]);
            gw.w.push_back(gw.w[i]);
        }
        
        std::reverse(gw.x.begin(),gw.x.end());
        
        gaussWeightsStorage.insert(std::pair<int, GaussWeight>(n, gw));
//        gaussWeightsStorage[n] = gw;

        x.resize(gw.x.size());
        w.resize(gw.x.size());
        
        x.assign(gw.x.begin(), gw.x.end());
        w.assign(gw.w.begin(), gw.w.end());
    } else {
        x.resize(it->second.x.size());
        w.resize(it->second.w.size());
        
        x.assign(it->second.x.begin(), it->second.x.end());
        w.assign(it->second.w.begin(), it->second.w.end());
    }
    pthread_mutex_unlock(&changeWeightsMutex);
};