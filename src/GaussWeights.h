//
//  GaussWeights.h
//  
//
//  Created by Andrew Nelson on 30/01/13.
//
//

#ifndef ____GaussWeights__
#define ____GaussWeights__

#include <iostream>
#include <vector>
#include <map>


class GaussWeight{
public:
    int n;
    std::vector<double> x;
    std::vector<double> w;
    GaussWeight();
    ~GaussWeight();
    
};


class GaussWeights{
        
public:
    GaussWeights();
    ~GaussWeights();
    void getGaussWeight(int n, std::vector<double> &x, std::vector<double> &w);
    
private:
    std::map<int,GaussWeight> gaussWeightsStorage;

};

#endif /* defined(____GaussWeights__) */
