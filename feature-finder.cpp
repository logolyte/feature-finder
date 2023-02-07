#include <iostream>

using namespace std;

#include <libxml2/libxml/parser.h>
#include <libxml2/libxml/tree.h>
#include <string.h>
#include <vector>
#include <fstream>


void print_tree(xmlNode* root, int generation) {
    for(int i=0; i<generation; i++) {
        cout<<' ';
    }
    cout<<root->name<<'\n';
    //if(generation>=11) return;
    xmlNode *child=root->children;
    if(child->next) {
        child=child->next;
        for(; child; child=child->next) {
            print_tree(child,generation+1);
        }
    }

}

vector<pair<unsigned int,unsigned int>> enhancer_indices;
vector<pair<unsigned int,unsigned int>> silencer_indices;

string read_sequence(string file_name) {
    string sequence;
    fstream filestream(file_name);
    string line;
    getline(filestream,line);
    while(getline(filestream,line)) {
        sequence.append(line);
    }
    filestream.close();
    return sequence;
}

void write_sequences(vector<pair<unsigned int,unsigned int>> index_vector, string sequence, string filename) {
    ofstream new_file(filename);
    pair<unsigned int,unsigned int> indices= *(index_vector.begin());
        string substring=sequence.substr(indices.first,indices.second-indices.first);
        new_file<<substring;
    for(vector<pair<unsigned int,unsigned int>>::iterator i=index_vector.begin()+1; i!=index_vector.end(); ++i) {
        pair<unsigned int,unsigned int> indices=*i;
        string substring=sequence.substr(indices.first,indices.second-indices.first);
        new_file<<'\n'<<substring;
    }
}

void get_indices(xmlNode* root) {
    vector<xmlNode*> node_list;
    vector<xmlNode*> next_list;
    node_list.push_back(root);
    while(node_list.size()>0) {
        for(auto i=node_list.begin(); i!=node_list.end(); ++i) {
            xmlNode *node=*i;
            if(node->content) {
                if(!strcmp((char*)node->content,"regulatory_class")) {
                    if(node->parent->next->next->children) {
                        string regulatory_type=(char*)node->parent->next->next->children->content;
                        if(!regulatory_type.compare("enhancer")) {
                            //cerr<<"enhancer, start: ";
                            xmlNode *location_node = node->parent->parent->parent->prev->prev->children->next->children->next->children->next;
                            int start = stoi((char*)location_node->children->next->children->content);
                            int stop =  stoi((char*)location_node->children->next->next->next->children->content);
                            //cerr<<start<<", stop: "<<stop<<'\n';
                            enhancer_indices.push_back(pair(start,stop));
                            //print_tree(location_node->children->next->next->next->next->next->children,100);
                            //cout<<location_node->children->next->next->next->next->next->children->next->name<<'\n';
                        } else if (!regulatory_type.compare("silencer")) {
                            //cerr<<"silencer, start: ";
                            xmlNode *location_node = node->parent->parent->parent->prev->prev->children->next->children->next->children->next;
                            int start = stoi((char*)location_node->children->next->children->content);
                            int stop =  stoi((char*)location_node->children->next->next->next->children->content);
                            //cerr<<start<<", stop: "<<stop<<'\n';
                            silencer_indices.push_back(pair(start,stop));
                        }
                    }
                }
            }
            for(xmlNode *child=node->children; child; child=child->next) {
                next_list.push_back(child);
            }
        }
        node_list.clear();
        for(auto i=next_list.begin(); i!=next_list.end(); ++i) {
            node_list.push_back(*i);
        }
        next_list.clear();
    }
}

void find_regulatory(xmlNode* root) {
    if(!strcmp((char *)root->name,"Gb-qual")) {
        if(root->next->content) {
            cout<<"found: "<<root->next->content<<'\n';
        }
    }
    xmlNode *child = root->children;
    for(; child; child=child->next) {
        find_regulatory(child);
    }

}

void analyze_file(string file_name){
    ifstream filestream;
    filestream.open(file_name);
    vector<string> enhancer_vector;
    string enhancer;
    int min_length=0;
    int max_length=0;
    while(getline(filestream,enhancer)){
        if(!min_length||enhancer.length()<min_length) min_length=enhancer.length();
        if(enhancer.length()>max_length) max_length=enhancer.length();
        enhancer_vector.push_back(enhancer);
    }
    cout<<"Number of entries: "<<enhancer_vector.size()<<'\n';
    cout<<"Minimum length: "<<min_length<<"\nMaximum length: "<<max_length<<'\n';
}

int main() {
    cout<<"Input XML file name:\n";
    string xml_name; //home/love/human-genome/chromosome-1.xml
    cin>>xml_name;
    xmlDoc *doc = xmlReadFile(xml_name.c_str(),NULL,0);
    xmlNode *root_element = xmlDocGetRootElement(doc);
    cout<<"Loading XML file succeeded\n";
    get_indices(root_element);
    cout<<"Getting enhancer and silencer indices succeeded\n";
    cout<<"Input fasta file name:\n";
    string fasta_name; //home/love/human-genome/chromosome-1.fasta
    cin>>fasta_name;
    string sequence=read_sequence(fasta_name);
    write_sequences(enhancer_indices,sequence,"enhancers.txt");
    cout<<"Wrote enhancers to the file enhancers.txt\n";
    write_sequences(silencer_indices,sequence,"silencers.txt");
    cout<<"Wrote silencers to the file silencers.txt\n";
    xmlFreeDoc(doc);
    xmlCleanupParser();
    cout<<"\nEnhancer summary: \n";
    analyze_file("enhancers.txt");
    cout<<"\nSilencer summary: \n";
    analyze_file("silencers.txt");
    cout << "Done" << endl;
    return 0;
}
