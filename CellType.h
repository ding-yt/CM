//
//  CellType.h
//  CancerModel_1
//
//  Created by Yuantong Ding on 11/18/13.
//  Copyright (c) 2013 Yuantong Ding. All rights reserved.
//

#ifndef __CancerModel_1__CellType__
#define __CancerModel_1__CellType__

#include <iostream>
#include <vector>

class Trait
{
protected:
    int name;
    double proliferation_time;
    double death_rate;
    double mutation_rate;
    double migration_rate;
    double fittness;
    
};

class Type : public Trait
{
public:
    Type();
    Type (int name_);
    void set_proliferationTime (double time);
    void set_deathRate (double rate);
    void set_mutationRate (double rate);
    void set_migrationRate (double rate);
    void set_fittness (double fit){fittness = fit;};
    
    int get_name() const { return name; }
    double get_proliferation_time() const { return proliferation_time; }
    double get_death_rate() const { return death_rate; }
    double get_mutation_rate() const { return mutation_rate; }
    double get_migration_rate() const { return migration_rate; }
    double get_fittness () const {return fittness;}
    
    void show();
    
};


class CellType
{
    std::vector<Type> all_types;
    int type_number;
    
public:
    CellType (int n_type);
    void set_typeNumber(int n_type){type_number = n_type;};
    Type& operator[] (const int nIndex);
    void show();
    double max_P();
};


#endif /* defined(__CancerModel_1__CellType__) */
