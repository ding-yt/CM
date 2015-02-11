//
//  Topology.cpp
//  NodeTester
//
//  Created by Yuantong Ding on 1/27/15.
//  Copyright (c) 2015 Yuantong Ding. All rights reserved.
//

#include "Topology.h"
#include <cassert>

Nodetree::Nodetree()
{
    _parent = NULL;
    _index = -1;
    _type = -1;
    _weight = 1;
    _name = "u";
    _sequence = "";
    
}

Nodetree::Nodetree(int index, int type, Nodetree * parent)
{
    _parent = parent;
    _index = index;
    _type = type;
    _weight = 1;
    _name = "u";
    _sequence = "";
}

Nodetree::~Nodetree(){
    
    for (int i=0; i<_children.size(); i++) {
        delete _children[i];
    }
    
}

void Nodetree::addName(int type, int x, int y){
    _name = "t";
    _name.append(std::to_string(type));
    _name.append("_");
    _name.append(std::to_string(x));
    _name.append("_");
    _name.append(std::to_string(y));
}

Topology::Topology():
_rng(ginkgo::RandomNumberGenerator::get_instance())
{
    _root = NULL;
    
}

Topology::~Topology(){
    destoryTree(_root);
}

void Topology::destoryTree(Nodetree * node){
    if (node != NULL){
        for (int i=0; i<node->_children.size(); i++) {
             delete node->_children[i];
        }
        delete node;
    }
}

Nodetree * Topology::searchNode(int index){
    //std::cout<<"root index:"<<_root->_index<<"\n";
    return searchNode(index, _root);
}

Nodetree * Topology::searchNode(int index, Nodetree * startNode){
    //std::cout<<"start Node index :"<<startNode->_index<<"\n";
    if (startNode==NULL) {
        return NULL;
    }
    if (startNode->_index==index){
        //std::cout<<"found\n";
        return startNode;
    }else{
        for (int i=0; i<startNode->_children.size(); i++) {
            if (startNode->_children[i]->_index == index) {
                return startNode->_children[i];
            }
        }
        return NULL;
    }
}

Nodetree * Topology::searchChildren(int index, Nodetree * parentNode){
    if (parentNode == NULL) {
        return NULL;
    }
    for (int i=0; i<parentNode->_children.size(); i++) {
        if (parentNode->_children[i]->_index == index) {
            return parentNode->_children[i];
        }
    }
    return NULL;

}


Nodetree * Topology::addNode(Nodetree * parentNode, int index, int type){    
    if (parentNode == NULL) {
        parentNode = new Nodetree(index, type, NULL);
        return parentNode;
    }
    Nodetree * children = searchChildren(index, parentNode);
    if (children == NULL) {
        Nodetree * newchild = new Nodetree(index,type,parentNode);
        parentNode->_children.push_back(newchild);
        return parentNode->_children[parentNode->_children.size()-1];
    }else{        
        return children;
    }
    
}

void Topology::addLineage(std::vector<CellIndexType> lineage){
    //std::cout <<lineage.size()<<"\n";
    if (lineage.size()<2) {
        return;
    }
    int endIndex = lineage.size()-1;
    //Nodetree * current = addNode(_root, lineage[endIndex-1], lineage[endIndex]);
    //_root = current;
    Nodetree * current = searchNode(lineage[endIndex-1]);
    for (int i=endIndex-2; i>=0; i-=2) {
        //std::cout<<"add cell:"<<lineage[i-1]<<" type: "<<lineage[i]<<"\n";
        current = addNode(current, lineage[i-1], lineage[i]);
//        if (i==1) {
//            current->addName(lineage[i], lineage[endIndex-2], lineage[endIndex-1]);
//        }
    }
}

void Topology::printNWK(){
    print(_root);
    std::cout<<"\n";
}

void Topology::print(Nodetree * node){
    if (node ==NULL) {
        std::cout<<"empty tree\n";
        return;
    }
    if (node->_children.size()==0) {
        //std::cout<<node->_index<<":"<<node->_weight;
        std::cout<<node->_index<<":"<<node->_weight;
        return;
    }
    std::cout<<"(";
    for (int i=0; i<node->_children.size(); i++) {
            print(node->_children[i]);
        if (i!=node->_children.size()-1) {
            std::cout<<",";
        }
        
//        }
    }
    std::cout<<")"<<node->_type<<":"<<node->_weight;
    
}

void Topology::addRoot(int index, int type){
    _root = new Nodetree(index,type,NULL);
}

void Topology::addChild(Nodetree *child){
    if (_root==NULL) {
        std::cout<<"no root!\n";
        return;
    }else{
        _root->_children.push_back(child);
    }
}

void Topology::compress(){
    compress(_root);
}

void Topology::compress(Nodetree * node){
    if (node->_children.size()==0) {
        return;
    }
    if (node->_children.size()==1 && node->_parent != NULL && node->_type == node->_children[0]->_type) {
        Nodetree * parent = node->_parent;
        Nodetree * child = node->_children[0];
        child->_weight += node->_weight;
        child->_parent = parent;
        for (int i=0; i<parent->_children.size(); i++) {
            if (parent->_children[i]->_index == node->_index) {
                parent->_children[i] = child;
            }
        }
        node->_parent = NULL;
        node->_children[0] = NULL;
        delete node;
        compress(child);
    }else{
        for (int i=0; i<node->_children.size(); i++) {
            compress(node->_children[i]);
        }
    }
}

void Topology::simulateSeq(std::vector<double> rate,int length){
    _root->_sequence = randomBase(length);
    int total_mutations = mutateSeq(_root, rate, 0);
    if (length>total_mutations) {
        chopSeq(_root,total_mutations);
    }
    
}

void Topology::chopSeq(Nodetree * node,int length){
    if (node->_children.size()==0) {
        node->_sequence = node->_sequence.substr(0,length);
        return;
    }
    
     node->_sequence = node->_sequence.substr(0,length);
    for (int i=0; i<node->_children.size(); i++) {
        chopSeq(node->_children[i], length);
    }
    
}

int Topology::mutateSeq(Nodetree * node, std::vector<double> rates, int count){
    
    if (node->_children.size()==0) {
        return count;
    }
    if (count>=node->_sequence.size()) {
        std::cout<<"warning: initial seq length too short!\n";
        return 0;
    }
    
    for (int i=0; i<node->_children.size(); i++) {
        node->_children[0]->_sequence = node->_sequence;
        float r = node->_children[i]->_weight * rates[node->_children[i]->_type];
        int m = _rng.poisson(r);
        for (int j=0; j<m; j++) {
            node->_children[i]->_sequence[count+j] = pointMutate(node->_children[i]->_sequence[count+j]);
        }
        count += m;
    }
    for (int i=0; i<node->_children.size(); i++) {
        mutateSeq(node->_children[i], rates, count);
    }
    
};

std::string Topology::randomBase(int n){
    long r;
    std::vector<std::string> bases{"A","T","C","G"};
    std::string s;
    for (int i=0;i<n; i++) {
        r = _rng.uniform_int(0, 3);
        s += bases[r];
    }
    return s;
}
char Topology::pointMutate(char s){
    std::vector<char> bases;
    if (s == 'A'){
        bases = {'T','C','G'};
    }else if (s=='T'){
        bases = {'A','C','G'};
    }else if (s=='C'){
        bases = {'A','T','G'};
    }else if (s=='G'){
        bases = {'A','C','T'};
    }else{
        std::cout<<"warning: unknown base "<<s<<"\n";
    }
    long r = _rng.uniform_int(0, 2);
    return bases[r];
    
}
