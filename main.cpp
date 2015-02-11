//
//  main.cpp
//  CM
//
//  Created by Yuantong Ding on 12/12/14.
//  Copyright (c) 2014 Yuantong Ding. All rights reserved.
//

#include <iostream>
#include "ParFile.h"
#include "Lattices.h"
#include "GenealogyNode.h"
#include "Topology.h"

int main(int argc, const char * argv[])
{
    
    std::string setting_file, output;
    std::map<std::string, double> parameters;
    ParFile p;
    int fail_count = 0;
    long random_seed;
    bool repeat_simulation = true;
    bool first_migration_event = true;
    map<double, vector<int>> migration_event;
    int migration_count = 0;
    int fileIndex = 0;
    
    clock_t start_time = clock();
    
    for (int i=1; i<argc; i++) {
        if (i+1 != argc) {
            if (strcmp(argv[i], "-i") == 0){
                setting_file = argv[i + 1];
            }else if (strcmp(argv[i], "-o") == 0){
                output = argv[i+1];
                //cout<<output<<"\n";
            }else if (strcmp(argv[i], "-s") == 0){
                random_seed = std::atoi(argv[i+1]);
            }
        }
        
    }
    
    output.erase(std::remove(output.begin(), output.end(), '\n'), output.end());
    std::cout<<"output is "<<output<<"\n";
    
    parameters = p.get_parameters(setting_file);
    p.show();
    GenealogyNode::set_allCellType(parameters);
    GenealogyNode::show_allCellType();
    
    double time_max = parameters["time_max"];
    int xStep = parameters["xStep"];
    int yStep = parameters["yStep"];
    double sample_time_interval = parameters["sample_time_interval"];
    int sample_size = (int)parameters["sample_size"];
    double time = 0;
    std::vector<double> pointMutation;
    pointMutation.push_back(parameters["type_0_point_mutation_rate"]);
    pointMutation.push_back(parameters["type_1_point_mutation_rate"]);
    pointMutation.push_back(parameters["type_2_point_mutation_rate"]);
    pointMutation.push_back(parameters["type_3_point_mutation_rate"]);
    
    GenealogyNode::set_counter(0);
    Lattices space(xStep,yStep);
    space.setRandomSeed(random_seed);
    space.getBase("/Users/dyt/Dropbox/cancerEvolution/script/CM/testBase.txt");
    
    while (repeat_simulation) {
        space.proliferate(xStep/2, yStep/2, 1);
        for (time=2; time<time_max; time++) {
            //Boundary condition check
            if (space.hitBoundary()) {
                repeat_simulation = false;
                cout<<"hit boundary! time:"<<time<<"\n";
                goto output;
            }
            
            //all normal cell check
            if (space.allNormalCell()) {
                repeat_simulation = true;
                fail_count ++;
                cout<<"faill attemp: allNormalCell\t"<<fail_count<<"\ttime:"<<time<<"\taliveCell:"<<space.aliveCellNumber()<<"\n";
                //string snapshot_fileName = output+to_string(fail_count)+".txt";
                //space.snapshot(snapshot_fileName, time);
                break;
            }else{
                repeat_simulation = false;
            }
            
            //pick a random cell to proliferate
            vector<int> loc = space.randomAliveCell();
            space.proliferate(loc[0], loc[1], time);
            
            //sampling
            if ((int)time%(int)sample_time_interval == 0 && (!repeat_simulation)) {
                // prepare output file
                std::string filename = output;
                std::ostringstream s;
                std::ostringstream intervel;
                intervel << sample_time_interval;
                std::string new_filename = filename + s.str()+".txt";
                fileIndex ++;
                s<<fileIndex;
                new_filename = filename + s.str()+"_i"+intervel.str()+".txt";
                               
                //sample file
                string sample_fileName = output+"_sample"+s.str()+".txt";

                space.snapshot(new_filename, time);
                space.sampleTumor(sample_size,sample_fileName,time);
            }
            
        }
        space.clear();
    }
    
output:
    // prepare output file
    std::string filename = output;
    cout<<"time"<<time<<"\n";
    string sample_fileName = output+"_sample"+"END.txt";
    string snapshot_fileName = output+"END.txt";
    string sampleSection_fileName = output +"_sectionEnd.txt";
    space.snapshot(snapshot_fileName, time);
    
    space.sampleTumor(sample_size,sample_fileName,time);
//    space.sampleSection(sample_size, sampleSection_fileName, time);
//    space.sampleLayer(sample_size, output, time);
    clock_t temp_time = clock();
    double sec = (temp_time-start_time)/CLOCKS_PER_SEC;
    int min = (int)sec/60;
    sec = sec - min*60;
    cout << "\nRunning time: " << min <<"min"<<sec<<"sec\n";
    
    
    return 0;
    
}

