//
//  Lattices.cpp
//  CancerModel_cleanup
//
//  Created by Yuantong Ding on 2/26/14.
//  Copyright (c) 2014 Yuantong Ding. All rights reserved.
//

#include "Lattices.h"
#include <fstream>
#include <sstream>
#include "randgen.hpp"
#include "Topology.h"

#define PI 3.141592653

using namespace std;


Lattices::Lattices(int xStep, int yStep, int initialCellNumber):
root_normal(new GenealogyNode()),
root_tumor(new GenealogyNode()),
_xStep(xStep),
_yStep(yStep),
_centerX(xStep/2),
_centerY(yStep/2),
_r(sqrt(initialCellNumber/PI)),
_initialCellNumber(initialCellNumber),
_emptyCellFittness(0),
_hitBoundary(false),
_migrationEvent(0),
_tumorXmin(50),
_tumorYmin(50),
_tumorXmax(50),
_tumorYmax(50),
rng_(ginkgo::RandomNumberGenerator::get_instance())
{
    double r2 = initialCellNumber/3.14;
    _aliveCellCount = 0;
    _normalCellCount = 0;
    for (int i=0; i<_xStep; i++) {
        vector<GenealogyNode*> g;
        cells.push_back(g);
        for (int j=0; j<_yStep; j++) {
            GenealogyNode * node = new GenealogyNode();
            cells[i].push_back(node);
            double ii = i-_xStep/2;
            double jj = j-_yStep/2;
            if (ii*ii+jj*jj<=r2) {
                cells[i][j]->set_parent(root_normal);
                _aliveCellCount ++;
                _normalCellCount ++;
            }else{
                cells[i][j]->die(0);
            }
            
        }
    }
    cells[_centerX][_centerY]->set_parent(root_tumor);
    cells[_centerX][_centerY]->set_type(1);
    _normalCellCount --;
    _treeTumor->setRoot(1, 1);
    _treeNormal->setRoot(0, 0);
    cout<<"root_normal: "<<root_normal->get_cell_index()<<"\n";
    cout<<"root_tumor: "<<root_tumor->get_cell_index()<<"\n";
}

void Lattices::setEmptyCellFit(double fit){
    _emptyCellFittness = fit;
}

Lattices::Lattices(int xStep, int yStep):
_xStep(xStep),
_yStep(yStep),
_centerX(xStep/2),
_centerY(yStep/2),
_emptyCellFittness(0),
_r(0),
_tumorXmin(50),
_tumorYmin(50),
_tumorXmax(50),
_tumorYmax(50),
_hitBoundary(false),
_migrationEvent(0),
_aliveCellCount(xStep*yStep),
_normalCellCount(xStep*yStep-1),
_initialCellNumber(xStep*yStep),
rng_(ginkgo::RandomNumberGenerator::get_instance())
{
    root_normal = new GenealogyNode();
    root_tumor = new GenealogyNode();
    for (int i=0; i<_xStep; i++) {
        vector<GenealogyNode*> g;
        cells.push_back(g);
        for (int j=0; j<_yStep; j++) {
            GenealogyNode * node = new GenealogyNode();
            cells[i].push_back(node);
            cells[i][j]->set_parent(root_normal);
        }
    }
    root_tumor->set_type(1);
    cells[_centerX][_centerY]->set_parent(root_tumor);
    cells[_centerX][_centerY]->set_type(1);
    cout<<"root_normal: "<<root_normal->get_cell_index()<<"\n";
    cout<<"root_tumor: "<<root_tumor->get_cell_index()<<"\n";
    cout<<"tumor at ("<<_centerX<<","<<_centerY<<"),type "<<cells[_centerX][_centerY]->get_type()<<"\n\n";
}

Lattices::~Lattices(){
    for (int i=0; i<_xStep; i++) {
        for (int j=0; j<_yStep; j++) {
            delete cells[i][j];
            cells[i][j]=NULL;
        }
    }
    root_normal = NULL;
    root_tumor = NULL;
}

void Lattices::clear(){
    for (int i=0; i<_xStep; i++) {
        for (int j=0; j<_yStep; j++) {
            delete cells[i][j];
            cells[i][j]=NULL;
        }
    }
    GenealogyNode::set_counter(0);
    root_normal = new GenealogyNode();
    root_tumor = new GenealogyNode();
//    root_normal->set_cell_index(0);
//    root_tumor->set_cell_index(1);
    for (int i=0; i<_xStep; i++) {
        for (int j=0; j<_yStep; j++) {
            cells[i][j]= new GenealogyNode();
            cells[i][j]->set_parent(root_normal);
        }
    }
    cells[_centerX][_centerY]->set_parent(root_tumor);
    cells[_centerX][_centerY]->set_type(1);
    _aliveCellCount = _xStep*_yStep;
    _normalCellCount = _aliveCellCount-1;
    _tumorXmin = 50;
    _tumorYmin = 50;
    _tumorXmax = 50;
    _tumorYmax = 50;
}

void Lattices::setRandomSeed(unsigned long seed){
    rng_.set_seed(seed);
}

void Lattices::printType(){
    std::cout<<"cell type:\n";
    for (int i=0; i<_xStep; i++) {
        for (int j=0; j<_yStep; j++) {
            cout<<cells[i][j]->get_type()<<"\t";
        }
        cout <<"\n";
    }
}

void Lattices::printCell(){
    for (int i=0; i<_xStep; i++) {
        for (int j=0; j<_yStep; j++) {
            cout<<cells[i][j]->get_cell_index()<<":"<<cells[i][j]->get_type()<<"\t";
        }
        cout <<"\n";
    }
}

GenealogyNode * Lattices::getCell(int x, int y){
    return cells[x][y];
}


/*
 * chose a random alive cell location from the lattices, return vector<int> loc, where loc[0] and loc[1]
 * are x and y coordinates respectively
 */
vector<int> Lattices::randomAliveCell(){
    vector<int> location;
    int rand = rng_.uniform_int(0, _xStep-1);
    location.push_back(rand);
    location.push_back(rng_.uniform_int(0, _yStep-1));
    while (!cells[location[0]][location[1]]->isAlive()) {
        location[0] = rng_.uniform_int(0, _xStep-1);
        location[1] = rng_.uniform_int(0, _yStep-1);
    }
    return location;
}

//////////////////??????????????????????????????
vector<int> Lattices::randomSampleAreaLimit(){
    vector<int> location;
    int original_r = sqrt(floor(_initialCellNumber/PI));
    int lowerX = _centerX-original_r;
    int upperX = _centerX+original_r;
    int lowerY = _centerY-original_r;
    int upperY = _centerY+original_r;
    int newR = _r;
    int rand = rng_.uniform_int(lowerX,upperX-1);
    
    location.push_back(rand);
    location.push_back(rng_.uniform_int(lowerY,upperY-1));
    double ii = location[0]-_centerX;
    double jj = location[1]-_centerY;
    double temp_r = ii*ii+jj*jj;
    
    //cout<<"random sample : "<<location[0]<<" "<<location[1]<<"original_r "<<original_r<<" _r "<<_r<<"\n";
    if (newR<original_r) {
        while (!cells[location[0]][location[1]]->isAlive() || temp_r>original_r*original_r) {
            location[0] = rng_.uniform_int(lowerX,upperX-1);
            location[1] = rng_.uniform_int(lowerY,upperY-1);
            ii = location[0]-_centerX;
            jj = location[1]-_centerY;
            temp_r = ii*ii+jj*jj;
        }
    }else{
        
        lowerX = _centerX-newR;
        upperX = _centerX+newR;
        lowerY = _centerY-newR;
        upperY = _centerY+newR;
        while (!cells[location[0]][location[1]]->isAlive() || temp_r>newR*newR) {
            location[0] = rng_.uniform_int(lowerX,upperX-1);
            location[1] = rng_.uniform_int(lowerY,upperY-1);
            ii = location[0]-_centerX;
            jj = location[1]-_centerY;
            temp_r = ii*ii+jj*jj;
        }
    }
    
    return location;
}

/*
 * print out and return cell lineage at a certain location, if cell is dead, print/return -1
 * lineage is a list of cell index of its parents and type
 * end with position of this cell
 * lineage: index of last child, type of last child, index of parents, type of parents,..., x, y
 */
vector<CellIndexType> Lattices::getLineage(int x,int y){
    vector<CellIndexType> indexList;
    if (!cells[x][y]->isAlive()) {
        std::cout<<"\t*\n";
        indexList.push_back(-1);
        return indexList;
    }else{
        GenealogyNode* parent = cells[x][y];
        while (parent->get_parent()) {
            //            std::cout<<"\t"<<parent->get_cell_index();
            indexList.push_back(parent->get_cell_index());
            indexList.push_back(parent->get_type());
            parent = parent->get_parent();
        }
        //        std::cout<<"\t"<<parent->get_cell_index()<<"\n";
        indexList.push_back(parent->get_cell_index());
        indexList.push_back(parent->get_type());
        indexList.push_back(x);
        indexList.push_back(y);
        return indexList;
    }
}


/*
 * return the location of a random neighbour (4 neighbors) of a cell at certain location
 */
vector<int> Lattices::randomNeighbour(int x, int y,double t){
    vector<int> neighbors = {x-1,y,x,y+1,x,y-1,x+1,y};
    vector<int> neighbor;
    for (int i=0; i<8; i+=2) {
        if (neighbors[i]>=0 && neighbors[i]<_xStep && neighbors[i+1]>=0 && neighbors[i+1]<_yStep) {
            neighbor.push_back(neighbors[i]);
            neighbor.push_back(neighbors[i+1]);
        }
    }
    vector<int> randomNeighbourLocation;
    if (neighbor.size()>0) {
        int rand = rng_.uniform_int(0,neighbor.size()/2-1);        
        randomNeighbourLocation.push_back(neighbor[rand*2]);
        randomNeighbourLocation.push_back(neighbor[rand*2+1]);
        //cout<<"Neighbour size: "<<neighbor.size()<<"\t random choice: "<<rand<<"\n";
    }
    

    return randomNeighbourLocation;
}


/*
 * return the location of a random neighbour(8 neighbors) of a cell at certain location
 */
vector<int> Lattices::randomNeighbour2(int x, int y,double t){
    vector<int> neighbors = {x-1,y,x,y+1,x,y-1,x+1,y,x+1,y+1,x+1,y-1,x-1,y+1,x-1,y-1};
    vector<int> neighbor;
    for (int i=0; i<8; i+=2) {
        if (neighbors[i]>=0 && neighbors[i]<_xStep && neighbors[i+1]>=0 && neighbors[i+1]<_yStep) {
            neighbor.push_back(neighbors[i]);
            neighbor.push_back(neighbors[i+1]);
        }
    }
    int rand = rng_.uniform_int(0,neighbor.size()/2-1);
    vector<int> randomNeighbourLocation;
    randomNeighbourLocation.push_back(neighbor[rand*2]);
    randomNeighbourLocation.push_back(neighbor[rand*2+1]);
    //    cout<<"Neighbour size: "<<neighbour_count<<"\t random choice: "<<rand<<"\n";
    return randomNeighbourLocation;
}


/*
 * decide a cell at certain location will die('D') or Migrate('M') or proliferate ('P')
 */
char Lattices::decideFate(int x,int y,double t){
    if (!cells[x][y]->isAlive()) {
        return 'D';
    }
    double rand = rng_.uniform_01();
    double d = cells[x][y]->get_death_rate();
    double m = cells[x][y]->get_death_rate()+cells[x][y]->get_migration_rate();
//    double sum_rate = cells[x][y]->get_death_rate()+cells[x][y]->get_migration_rate()+cells[x][y]->get_proliferation_time()/cells[x][y]->get_maxP();
//    double d = cells[x][y]->get_death_rate()/sum_rate;
//    double p = (d+cells[x][y]->get_proliferation_time()/cells[x][y]->get_maxP())/sum_rate;
    //cout<<"cell at ("<<x<<", "<<y<<"), type "<<cells[x][y]->get_type()<<", d:"<<cells[x][y]->get_death_rate()<<", m"<<cells[x][y]->get_migration_rate()<<", p:"<<cells[x][y]->get_proliferation_time()/cells[x][y]->get_maxP()<<"\n";
    //cout<<sum_rate<<"\t"<<d<<"\t"<<p<<"\n";


    
    if (rand<d) {
        //cout<<rand<<" D\n";
        return 'D';
    }else if(rand<m){
        //cout<<rand<<" P\n";
        return 'M';
    }else{
        //cout<<d<<"\t"<<p<<"\tdie\t"<<rand<<"\n";
        //cout<<rand<<" M\n";
        return 'P';
    }
    
}

/*
 * For a cell at location (x,y), decide whether it will migrate
 */
bool Lattices::migrate(int x, int y, double t){
    if (!cells[x][y]->isAlive()) {
        return false;
    }
    double rand = rng_.uniform_01();
    //    cout<<"Migration rate: "<<cells.get(x, y)->get_migration_rate()<<"\t random: "<<rand<<"\n";//
    if (rand<cells[x][y]->get_migration_rate()) {
        cells[x][y]->die(t);
        _aliveCellCount --;
        if (cells[x][y]->get_type()==0) {
            _normalCellCount --;
        }
        return true;
    }else{
        return false;
    }
}

/*
 * For a cell at location (x,y) go through proliferation, decide whether it will mutate
 */
bool Lattices::mutate(int x, int y, double t){
    if (!cells[x][y]->isAlive()) {
        return false;
    }
    double rand = rng_.uniform_01();
    //    cout<<"Mutation rate: "<<cells.get(x, y)->get_mutation_rate()<<"\t random: "<<rand<<"\n";//
    if (rand<cells[x][y]->get_mutation_rate()) {
        return true;
    }else{
        return false;
    }
}

bool Lattices::replace(double fit, int n, int m){
    
    if (cells[n][m]->isAlive()) {
        double f = fit/(fit+cells[n][m]->get_fittness());
        double rand = rng_.uniform_01();
        //cout<<"fitness "<<fit<<" replace "<<cells[n][m]->get_fittness()<<", ratio "<<f<<", rand "<<rand<<"\n";
        if (rand<f) {
            //cout<<"\t sucess\n";
            return true;
        }else{
            return false;
        }
    }else{
        return true;
    }
}


/*
 * For a cell at location (x,y), return its dead neighbour location as x-y key-value pair
 * if there's no dead neighbour, map size is 0
 */
vector<int> Lattices::emptyNeighbour(int x,int y,double t){
    vector<int> emptyNeighbours;
    vector<int> neighbors = {x-1,y,x,y+1,x,y-1,x+1,y,x+1,y+1,x+1,y-1,x-1,y+1,x-1,y-1};
    vector<int> neighbor;
    for (int i=0; i<4; i+=2) {
        if (neighbors[i]>0 && neighbors[i]<_xStep && neighbors[i+1]>0 && neighbors[i+1]<_yStep && !cells[neighbor[i]][neighbor[i+1]]->isAlive()) {
            neighbor.push_back(neighbors[i]);
            neighbor.push_back(neighbors[i+1]);
        }
    }
    return emptyNeighbours;
}

/*
 * For a cell at location (x,y), if go through proliferation, decide whether its 2 offsprings mutate
 * The first offspring will occupy its random dead neighbour's space if available, or compete with a random alive neighbour
 * The second offspring will replace the parent's space
 */
void Lattices::proliferate(int x,int y, double t){
    char stage = decideFate(x, y, t);
    //cout<<"cell at ("<<x<<","<<y<<"), type "<<cells[x][y]->get_type()<<" stage "<<stage<<"\n";
    if (stage == 'P') {        
        double p = cells[x][y]->get_proliferation_time()/cells[x][y]->get_maxP();
        double rand = rng_.uniform_01();
        //cout<<"\tp "<<p<<" rand:"<<rand<<"\n";
        if (rand<p) {            
            // offspring mutate?
            GenealogyNode * temp_offsprint;
            if (mutate(x, y, t)) {
                temp_offsprint = new GenealogyNode(cells[x][y]->get_type()+1,t);
            }else{
                temp_offsprint = new GenealogyNode(cells[x][y]->get_type(),t);
            }
            
            vector<int> neighbor = randomNeighbour(x, y, t);
           
            
            if (neighbor.size()>0) {
         
            int n_x = neighbor[0];
            int n_y = neighbor[1];
            
             //cout<<"\t p success, offspring type "<<temp_offsprint->get_type()<<"\n";
                //cout<<"\t replace cell at ("<<n_x<<","<<n_y<<"), type "<<cells[n_x][n_y]->get_type()<<"\n";
                
            if (replace(temp_offsprint->get_fittness(), neighbor[0], neighbor[1])) {
                temp_offsprint->set_parent(cells[x][y]);
                if ( n_x==0 || n_x==_xStep-1 || n_y==0 || n_y==_yStep-1) {
                    if (temp_offsprint->get_type()>0) {
                        _hitBoundary = true;
                    }
                }
                if (cells[n_x][n_y]->isAlive()) {
                    if (cells[n_x][n_y]->get_type()==0 && temp_offsprint->get_type()>0) {
                        _normalCellCount --;
                    }
                    if (cells[n_x][n_y]->get_type()>0 && temp_offsprint->get_type()==0) {
                        _normalCellCount ++;
                    }
                    cells[n_x][n_y]->die(t);
                }else{
                    if (temp_offsprint->get_type()==0) {
                        _normalCellCount ++;
                    }
                    _aliveCellCount ++;
                }
                
                cells[n_x][n_y] = temp_offsprint;
                if (cells[n_x][n_y]->get_type()>0) {
                    if (n_x<_tumorXmin) {
                        _tumorXmin = n_x;
                    }else if (n_x>_tumorXmax){
                        _tumorXmax = n_x;
                    }
                    if (n_y<_tumorYmin) {
                        _tumorYmin = n_y;
                    }else if (n_y > _tumorYmax){
                        _tumorYmax = n_y;
                    }
                }
                
            }else{
                delete temp_offsprint;
            }
            }
            
        
      
        // self mutate?
        //        cout<<"self mutate?\n";//
        GenealogyNode * tempself;
        if (mutate(x, y, t)) {
            if (cells[x][y]->get_type()==0) {
                _normalCellCount --;
            }
            tempself = new GenealogyNode(cells[x][y]->get_type()+1,t);
          
        }else{
            tempself = new GenealogyNode(cells[x][y]->get_type(),t);
        }
        tempself->set_parent(cells[x][y]);
        cells[x][y]= tempself;
            
        }
        
    }else{  //'D' or 'M'
        //cout<<"die\n";
        if(stage=='M'){
            _migrationEvent ++;
        }
        if (cells[x][y]->isAlive()) {
            if (cells[x][y]->get_type()==0) {
                _normalCellCount --;
            }
            //cells[x][y]->die(t);
            _aliveCellCount --;
        }
    }
}


bool Lattices::allNormalCell(){
    //    cout << "aliveCellCount: "<<_aliveCellCount<<" normalCellCount :"<<_normalCellCount<<"\n";
    if (_aliveCellCount == _normalCellCount) {
        return true;
    }else{
        return false;
    }
}

void Lattices::sampling(int sampleSize, string filename,double time){
    vector<int> locationX;
    vector<int> locationY;
    cout<<"samplefile: "<<filename<<"\n";
    for (int i=0; i<sampleSize; i++) {
        int x = rng_.uniform_int(0, _xStep-1);
        int y = rng_.uniform_int(0, _yStep-1);
        if (!cells[x][y]->isAlive()) {
            i--;
        }else{
            locationX.push_back(x);
            locationY.push_back(y);
        }
    }
    ofstream sample_file;
    sample_file.open(filename.c_str());
    
    int normalCount =0;
    int type1count = 0;
    int type2count = 0;
    int type3count = 0;
    for (int i=0; i<_xStep; i++) {
        for (int j=0; j<_yStep; j++) {
            if (cells[i][j]->get_type()==0) {
                normalCount++;
            }else if(cells[i][j]->get_type()==1){
                type1count++;
            }else if(cells[i][j]->get_type()==2){
                type2count++;
            }else if(cells[i][j]->get_type()==3){
                type3count++;
            }
        }
    }
    sample_file<<"sample at time "<<time<<"\n";
    sample_file<<"type_0_count\t"<<normalCount<<"\n";
    sample_file<<"type_1_count\t"<<type1count<<"\n";
    sample_file<<"type_2_count\t"<<type2count<<"\n";
    sample_file<<"type_3_count\t"<<type3count<<"\n";
    
    for (int s=0; s<sampleSize; s++) {
        int sampleCount = 0;
        int x = locationX[s];
        int y = locationY[s];
        vector<CellIndexType> lineage = getLineage(x, y);
//        if (lineage[lineage.size()-4]==0) {
//            _treeNormal->addLineage(lineage);
//        }else if(lineage[lineage.size()-4]==1){
//            _treeNormal->addLineage(lineage);
//        }
        sample_file<<"sample"<<s<<"\n";
        for (int i=0; i<lineage.size(); i++) {
            sample_file<<lineage[i]<<"\t";
        }
        sample_file<<"\n";
        
    }
    sample_file.close();
    
//    _treeNormal->compression();
//    _treeTumor->compression();
//    
//    _treeTumor->printTree();
    
    
    
}

void Lattices::sampleTumor(int sampleSize, string filename,double time){
    vector<int> locationX;
    vector<int> locationY;
    cout<<"samplefile: "<<filename<<"\n";
    for (int i=0; i<sampleSize; i++) {
        int x = rng_.uniform_int(_tumorXmin, _tumorXmax);
        int y = rng_.uniform_int(_tumorYmin, _tumorYmax);
        if (!cells[x][y]->isAlive()) {
            i--;
        }else{
            locationX.push_back(x);
            locationY.push_back(y);
        }
    }
    ofstream sample_file;
    sample_file.open(filename.c_str());
    
    int normalCount =0;
    int type1count = 0;
    int type2count = 0;
    int type3count = 0;
    for (int i=0; i<_xStep; i++) {
        for (int j=0; j<_yStep; j++) {
            if (cells[i][j]->get_type()==0) {
                normalCount++;
            }else if(cells[i][j]->get_type()==1){
                type1count++;
            }else if(cells[i][j]->get_type()==2){
                type2count++;
            }else if(cells[i][j]->get_type()==3){
                type3count++;
            }
        }
    }
    sample_file<<"sample at time "<<time<<"\n";
    sample_file<<"type_0_count\t"<<normalCount<<"\n";
    sample_file<<"type_1_count\t"<<type1count<<"\n";
    sample_file<<"type_2_count\t"<<type2count<<"\n";
    sample_file<<"type_3_count\t"<<type3count<<"\n";
    
    for (int s=0; s<sampleSize; s++) {
        int sampleCount = 0;
        int x = locationX[s];
        int y = locationY[s];
        vector<CellIndexType> lineage = getLineage(x, y);
//        if (lineage[lineage.size()-4]==0) {
//            _treeNormal->addLineage(lineage);
//        }else if(lineage[lineage.size()-4]==1){
//            _treeNormal->addLineage(lineage);
//        }
        sample_file<<"sample"<<s<<"\n";
        for (int i=0; i<lineage.size(); i++) {
            sample_file<<lineage[i]<<"\t";
        }
        sample_file<<"\n";
        
    }
    sample_file.close();
    
    nwk(sampleSize, locationX, locationY);
    
//    _treeNormal->compression();
//    _treeTumor->compression();
//    
//    _treeTumor->printTree();
    
}

void Lattices::snapshot(string filename,double t){
    cout<<"all cell: "<<filename<<"\n";
    ofstream snapshot_file;
    snapshot_file.open(filename.c_str());
    snapshot_file<<"cell at time "<<t<<"\n";
    snapshot_file<<"migration event count "<<_migrationEvent<<"\n";
    
    int normalCount =0;
    int type1count = 0;
    int type2count = 0;
    int type3count = 0;
    for (int i=0; i<_xStep; i++) {
        for (int j=0; j<_yStep; j++) {
            if (cells[i][j]->get_type()==0) {
                normalCount++;
            }else if(cells[i][j]->get_type()==1){
                type1count++;
            }else if(cells[i][j]->get_type()==2){
                type2count++;
            }else if(cells[i][j]->get_type()==3){
                type3count++;
            }
        }
    }
    snapshot_file<<"type_0_count\t"<<normalCount<<"\n";
    snapshot_file<<"type_1_count\t"<<type1count<<"\n";
    snapshot_file<<"type_2_count\t"<<type2count<<"\n";
    snapshot_file<<"type_3_count\t"<<type3count<<"\n";
    for (int i=0; i<_xStep; i++) {
        for (int j=0; j<_yStep; j++) {
            snapshot_file<<i<<"\t"<<j<<"\t";
            if (cells[i][j]->isAlive()) {
                snapshot_file<<cells[i][j]->get_type()<<"\n";
            }else{
                snapshot_file<<"-1\n";
            }
        }
    }
    snapshot_file.close();
    
//    string sfilename = filename;
//    sfilename.append(".symbol");
//    ofstream pic;
//    pic.open(sfilename.c_str());
//    for (int i=0; i<_xStep; i++) {
//        for (int j=0; j<_yStep; j++) {
//            
//            if (cells[i][j]->isAlive()) {
//                pic<<cells[i][j]->get_type();
//            }else{
//                pic<<".";
//            }
//            
//        }
//        pic<<"\n";
//    }
//    pic.close();
    
    
    string spfilename = filename.append(".symbol");
//    string spfilename = filename;
//    spfilename.append(".symbol");
    ofstream sp;
    sp.open(spfilename.c_str());
    for (int i=0; i<_xStep; i++) {
        for (int j=0; j<_yStep; j++) {
            if (i<601 && i >300 && j>300 && j<601 ) {
                switch (cells[i][j]->get_type()) {
                    case 1:
                        sp<<"A";
                        break;
                    case 2:
                        sp<<"B";
                        break;
                    case 3:
                        sp<<"C";
                        break;
                    default:
                        sp<<"N";
                        break;
                }
            
            }else{
                if (cells[i][j]->isAlive()) {
                    sp<<cells[i][j]->get_type();
                }else{
                    sp<<".";
                }
            }
        }
        sp<<"\n";
    }
    sp.close();
    
//    string lfilename = filename.append(".allLineage.txt");
//    //string lfilename = filename;
//    //lfilename.append(".allLineage.txt");
//    ofstream lineage;
//    lineage.open(lfilename.c_str());
//    lineage<<"snapshot at time "<<t<<"\n";
//    for (int i=301; i<601; i++) {
//        for (int j=301; j<601; j++) {
//            if (cells[i][j]->get_type()==1) {
//                vector<CellIndexType> lin = getLineage(i, j);
//                lineage<<i<<"\t"<<j<<"\t"<<cells[i][j]->get_type()<<"\n";
//                for (int i=0; i<lin.size(); i++) {
//                    lineage<<lin[i]<<"\t";
//                }
//                lineage<<"\n";
//            }
//
//        }
//        lineage<<"\n";
//    }
//    lineage.close();
}

void Lattices::sampleSection(int sampleSize, string filename,double time){
    vector<int> locationX;
    vector<int> locationY;
    cout<<"samplefile: "<<filename<<"\n";
    if (_tumorXmax-_tumorXmin>_tumorYmax-_tumorYmin) {
        int section = (_tumorXmax-_tumorXmin)/5;
        for (int i=0; i<5; i++) {
            for (int j=0; j<sampleSize/5; j++) {
                int x = rng_.uniform_int(_tumorXmin+i*section, _tumorXmin+(i+1)*section);
                int y = rng_.uniform_int(_tumorYmin, _tumorYmax);
                if (!cells[x][y]->isAlive()) {
                    i--;
                }else{
                    locationX.push_back(x);
                    locationY.push_back(y);
                }
            }
        }
    }else{
        int section = (_tumorYmax-_tumorYmin)/5;
        for (int i=0; i<5; i++) {
            for (int j=0; j<sampleSize/5; j++) {
                int x = rng_.uniform_int(_tumorXmin, _tumorXmax);
                int y = rng_.uniform_int(_tumorYmin+i*section, _tumorYmin+(i+1)*section);
                if (!cells[x][y]->isAlive()) {
                    i--;
                }else{
                    locationX.push_back(x);
                    locationY.push_back(y);
                }
            }
        }
    }
    
    ofstream sample_file;
    sample_file.open(filename.c_str());
    
    int normalCount =0;
    int type1count = 0;
    int type2count = 0;
    int type3count = 0;
    for (int i=0; i<_xStep; i++) {
        for (int j=0; j<_yStep; j++) {
            if (cells[i][j]->get_type()==0) {
                normalCount++;
            }else if(cells[i][j]->get_type()==1){
                type1count++;
            }else if(cells[i][j]->get_type()==2){
                type2count++;
            }else if(cells[i][j]->get_type()==3){
                type3count++;
            }
        }
    }
    sample_file<<"sample at time "<<time<<"\n";
    sample_file<<"type_0_count\t"<<normalCount<<"\n";
    sample_file<<"type_1_count\t"<<type1count<<"\n";
    sample_file<<"type_2_count\t"<<type2count<<"\n";
    sample_file<<"type_3_count\t"<<type3count<<"\n";
    
    for (int s=0; s<sampleSize; s++) {
        int sampleCount = 0;
        int x = locationX[s];
        int y = locationY[s];
        vector<CellIndexType> lineage = getLineage(x, y);
        sample_file<<"sample"<<s<<"\n";
        for (int i=0; i<lineage.size(); i++) {
            sample_file<<lineage[i]<<"\t";
        }
        sample_file<<"\n";
        
    }
    sample_file.close();
    
}

void Lattices::sampleLayer(int sampleSize, string filename,double time){
    vector<int> locationX_core;
    vector<int> locationY_core;
    vector<int> locationX_boundary;
    vector<int> locationY_boundary;
    
    int r;
    int x_trans = _tumorXmin+(_tumorXmax-_tumorXmin)/2;
    int y_trans = _tumorYmin+(_tumorYmax-_tumorYmin)/2;
    if (_tumorXmax-_tumorXmin>_tumorYmax-_tumorYmin) {
        r = (_tumorYmax-_tumorYmin)/2;
    }else{
        r = (_tumorXmax-_tumorXmin)/2;
    }
    int r2 = r*r/2;
    r = (int)r/sqrt(r);
    
    for (int i=0; i<sampleSize; i++) {
        int x = rng_.uniform_int(x_trans-r,x_trans+r);
        int y = rng_.uniform_int(y_trans-r,y_trans+r);

        if (!cells[x][y]->isAlive()) {
            i--;
        }else{
            int tx = x-x_trans;
            int ty = y-y_trans;
            if (tx*tx+ty*ty<=r2) {
                locationX_core.push_back(x);
                locationY_core.push_back(y);
            }
        }
    }
    
    for (int i=0; i<sampleSize; i++) {
        int x = rng_.uniform_int(_tumorXmin, _tumorXmax);
        int y = rng_.uniform_int(_tumorYmin, _tumorYmax);
        if (!cells[x][y]->isAlive()) {
            i--;
        }else{
            int tx = x-x_trans;
            int ty = y-y_trans;
            if (tx*tx+ty*ty>r2) {
                locationX_boundary.push_back(x);
                locationY_boundary.push_back(y);
            }else{
                i--;
            }
        }
    }

    
    
    
    ofstream sample_file;
    string coreSample = filename;
    coreSample.append("_coreEnd.txt");
    sample_file.open(coreSample.c_str());
    cout<<"samplefile: "<<coreSample<<"\n";
    
    int normalCount =0;
    int type1count = 0;
    int type2count = 0;
    int type3count = 0;
    for (int i=0; i<_xStep; i++) {
        for (int j=0; j<_yStep; j++) {
            if (cells[i][j]->get_type()==0) {
                normalCount++;
            }else if(cells[i][j]->get_type()==1){
                type1count++;
            }else if(cells[i][j]->get_type()==2){
                type2count++;
            }else if(cells[i][j]->get_type()==3){
                type3count++;
            }
        }
    }
    sample_file<<"sample at time "<<time<<"\n";
    sample_file<<"type_0_count\t"<<normalCount<<"\n";
    sample_file<<"type_1_count\t"<<type1count<<"\n";
    sample_file<<"type_2_count\t"<<type2count<<"\n";
    sample_file<<"type_3_count\t"<<type3count<<"\n";
    
    for (int s=0; s<sampleSize; s++) {
        int sampleCount = 0;
        int x = locationX_core[s];
        int y = locationY_core[s];
        vector<CellIndexType> lineage = getLineage(x, y);
        sample_file<<"sample"<<s<<"\n";
        for (int i=0; i<lineage.size(); i++) {
            sample_file<<lineage[i]<<"\t";
        }
        sample_file<<"\n";
        
    }
    sample_file.close();
    
    ofstream sample_file_boundary;
    string boundarySample = filename;
    boundarySample.append("_boundaryEnd.txt");
    sample_file_boundary.open(boundarySample.c_str());
    cout<<"samplefile: "<<boundarySample<<"\n";

    sample_file_boundary<<"sample at time "<<time<<"\n";
    sample_file_boundary<<"type_0_count\t"<<normalCount<<"\n";
    sample_file_boundary<<"type_1_count\t"<<type1count<<"\n";
    sample_file_boundary<<"type_2_count\t"<<type2count<<"\n";
    sample_file_boundary<<"type_3_count\t"<<type3count<<"\n";
    
    for (int s=0; s<sampleSize; s++) {
        int sampleCount = 0;
        int x = locationX_boundary[s];
        int y = locationY_boundary[s];
        vector<CellIndexType> lineage = getLineage(x, y);
        sample_file_boundary<<"sample"<<s<<"\n";
        for (int i=0; i<lineage.size(); i++) {
            sample_file_boundary<<lineage[i]<<"\t";
        }
        sample_file_boundary<<"\n";
        
    }
    sample_file_boundary.close();
    
    
}


void Lattices::nwk(int sampleSize,vector<int> locationX,vector<int> locationY){
    Topology * t = new Topology();
    Topology * tT = new Topology();
    Topology * tAll = new Topology();
    t->addRoot(0, 0);
    tT->addRoot(1, 1);
    tAll->addRoot(-1, -1);
    
    for (int s=0; s<sampleSize; s++) {
        int sampleCount = 0;
        int x = locationX[s];
        int y = locationY[s];
        vector<CellIndexType> lineage = getLineage(x, y);
        lineage.pop_back();
        lineage.pop_back();
        if (lineage[lineage.size()-2]==0) {
            t->addLineage(lineage);
        }else if (lineage[lineage.size()-2]==1){
            tT->addLineage(lineage);
        }
    }
    
    
    t->compress();
    tT->compress();
    
    tAll->addChild(t->getRoot());
    tAll->addChild(tT->getRoot());
    
//    tAll->printNWK();
//     
//    std::cout<<"Normal tree:\n";
//    t->printNWK();
    std::cout<<"\nTumor tree:\n";
    tT->printNWK();
}

//void Lattices::simulateSeq(std::vector<double> rate, Topology * t, int size){
//    t->setSeq(randomBase(size));
//    
//}


void Lattices::getBase(std::string filename){
    
    std::vector<std::vector<int>> base;
    std::ifstream file;
    std::string line;
    file.open(filename.c_str());
    if (file.is_open()) {
        
        while (std::getline(file, line)) {
            std::vector<int> baseLineage;
            std::stringstream ss(line);
            int n;
            if (ss>>n){
            while (ss >> n){
                baseLineage.push_back(-n);
            }
            base.push_back(baseLineage);
            }
        }
        
        file.close();
    }else{
        std::cout << "Unable to open file "<<filename<<"\n";
    }
    
    std::cout <<"lineage base\t"<<base.size()<<"\n";

    for (int i=0; i<base[0].size(); i++) {
        std::cout <<base[0][i]<<"\t";
    }
     std::cout <<"\n";
}

