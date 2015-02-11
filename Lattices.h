//
//  Lattices.h
//  CancerModel_cleanup
//
//  Created by Yuantong Ding on 2/26/14.
//  Copyright (c) 2014 Yuantong Ding. All rights reserved.
//

#ifndef __CancerModel_cleanup__Lattices__
#define __CancerModel_cleanup__Lattices__

#include <iostream>
#include <vector>
#include <map>
#include "GenealogyNode.h"
#include "randgen.hpp"
#include "Tree.h"
#include <math.h>
#include "Topology.h"


using namespace std;

class Lattices{
    std::vector<std::vector<GenealogyNode*>> cells;
    int _xStep;
    int _yStep;
    GenealogyNode * root_normal;
    GenealogyNode * root_tumor;
    ginkgo::RandomNumberGenerator& rng_;
    int _aliveCellCount;
    int _normalCellCount;
    double _emptyCellFittness;
    
    double _r;
    int _centerX;
    int _centerY;
    int _initialCellNumber;
    bool _hitBoundary;
    int _migrationEvent;
    int _tumorXmin;
    int _tumorYmin;
    int _tumorXmax;
    int _tumorYmax;
    Tree * _treeTumor;
    Tree * _treeNormal;
    
public:
    Lattices(int xStep,int yStep, int initialCellNumber);
    
    Lattices(int xStep, int yStep);
    
    ~Lattices();
    
    void setRandomSeed(unsigned long seed);
    
    void printType();
    
    void printCell();
    
    void setEmptyCellFit(double fit);
        
    char decideFate(int x,int y,double t);
    
    bool migrate(int x,int y, double t);
    
    bool mutate(int x,int y,double t);
    
    bool replace(double fit, int n, int m);
    
    vector<int> emptyNeighbour(int x,int y,double t);
    
    vector<int> randomNeighbour(int x,int y,double t);
    
    vector<int> randomNeighbour2(int x,int y,double t);
    
    void proliferate(int x,int y, double t);
    
    vector<int> randomAliveCell();
    
    vector<int> randomSampleAreaLimit();
    
    vector<CellIndexType> getLineage(int x,int y);
    
    GenealogyNode * getCell(int x,int y);
        
    int normalCellNumber(){return _normalCellCount;};
    
    int aliveCellNumber(){return _aliveCellCount;};
    
    bool allNormalCell();
    
    void updateOxygen();
    
    void clear();
    
    int getCenterX(){return _centerX;};
    
    int getCenterY(){return _centerY;};
    
    double getSampleR(){return _r;};
    
    bool hitBoundary(){return _hitBoundary;};
    
    void sampling(int sampleSize, string filename,double time);
    
    void sampleTumor(int sampleSize, string filename,double time);
    
    void sampleSection(int sampleSize, string filename,double time);
    
    void sampleLayer(int sampleSize, string filename,double time);
    
    void snapshot(string filename,double t);
    
    void nwk(int sampleSize,vector<int> X,vector<int> Y);
    
    void getBase(string filename);
    
//    void simulateSeq(std::vector<double> rate,Topology * t,int size);
    
//private:
//    string randomBase(int n);
//    string pointMutate(string s);

};


#endif /* defined(__CancerModel_cleanup__Lattices__) */
