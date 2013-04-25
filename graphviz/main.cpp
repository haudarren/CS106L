/******************************************************
 * File: main.cpp
 * Author: Darren Hau
 * -----
 * This program implements a Fruchterman-Reingold force-directed layout of a graph.  The program
 * prompts the user for a file, run time, and filter value.  The filter value takes in a specified
 * number of delta_x and delta_y to be averaged in a final delta_avg, simulating inertia.  The program
 * also uses random perturbations to avoid local optimizations.
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <deque>
#include "GraphVisualizer.h"
#include "SimpleGraph.h"
using namespace std;

/* Constants */
const double PI = 3.14159265358979323;
const double K_REPEL = 0.001;
const double K_ATTRACT = 0.001;

/* Function Prototypes */
void Welcome();
string promptUserForFile(string prompt);
double promptUserForTime(string prompt);
bool promptUserForLoop(string prompt);
int promptUserForFilter(string prompt);
string toUpperCase(string str);
void loadFileToGraph(SimpleGraph & graph, ifstream & stream);
void initNodePositions(SimpleGraph & graph, int numNodes);
void runFruchtermanReingold(SimpleGraph & graph, double requestedTime, int filterNum);
void computeRepelForces(SimpleGraph & graph, vector<double> & dx_tot, vector<double> & dy_tot, int node0, int node1);
void computeAttractForces(SimpleGraph & graph, vector<double> & dx_tot, vector<double> & dy_tot, int edgeNum);
void filterVelocities(vector<double> & dx_avg, vector<double> & dy_avg, vector< deque<double> > & dx_filter, vector< deque<double> > & dy_filter, vector<double> & dx_tot, vector<double> & dy_tot, int numNodes, int filterNum);

/* Main */
int main() {
	Welcome();
    ifstream infile;
    
    while (true) {
        string file = promptUserForFile("Enter file name: ");
        infile.open(file.c_str());
        
        double requestedTime = promptUserForTime("Enter time: ");
        int filterNum = promptUserForFilter("Enter a filter value between 1 and 500: ");

        SimpleGraph graph;
        loadFileToGraph(graph, infile);
        InitGraphVisualizer();
        runFruchtermanReingold(graph, requestedTime, filterNum);
        infile.close();
        
        bool loop = promptUserForLoop("Would you like to display another graph?");
        if (!loop) break;
    }
	return 0;
}

/* Function: runFruchtermanReingold
 * Usage: runFruchtermanReingold(graph, requestedTime, filterNum)
 * -----
 * Lays out a graph according to the Fruchterman-Reingold force-directed algorithm.
 */
void runFruchtermanReingold(SimpleGraph & graph, double requestedTime, int filterNum) {
    time_t startTime = time(NULL);
    int numNodes = graph.nodes.size();
    
    // This filter discards the oldest values and pushes in the newest value.
    vector< deque<double> > dx_filter(numNodes, deque<double>(filterNum));
    vector< deque<double> > dy_filter(numNodes, deque<double>(filterNum));
    
    while (true) {
        DrawGraph(graph);
        
        vector<double> dx_tot(numNodes);
        vector<double> dy_tot(numNodes);
        
        for (int i = 0; i < numNodes; i++) {            // compute repellent forces
            for (int j = i+1; j < numNodes; j++) {
                computeRepelForces(graph, dx_tot, dy_tot, i, j);
            }
        }
        
        for (int i = 0; i < graph.edges.size(); i++) {  // compute attractive forces
            computeAttractForces(graph, dx_tot, dy_tot, i);
        }
        
        vector<double> dx_avg(numNodes);
        vector<double> dy_avg(numNodes);
        filterVelocities(dx_avg, dy_avg, dx_filter, dy_filter, dx_tot, dy_tot, numNodes, filterNum);
        
        for (int i = 0; i < numNodes; i++) {            // update node positions
            double randomShift = (double)(rand() % 10) / RAND_MAX;
            graph.nodes[i].x += dx_avg[i] + randomShift;
            graph.nodes[i].y += dy_avg[i] - randomShift;
        }
        
        double elapsedTime = difftime(time(NULL), startTime);
        if (elapsedTime >= requestedTime) break;
    }
}

/* Function: filterVelocities
 * Usage: filterVelocities(dx_avg, dy_avg, dx_filter, dy_filter, dx_tot, dy_tot, numNodes, filterNum)
 * -----
 * Averages the past "filterNum" number of delta_x and delta_y to simulate inertia.
 */
void filterVelocities(vector<double> & dx_avg, vector<double> & dy_avg, vector< deque<double> > & dx_filter, vector< deque<double> > & dy_filter, vector<double> & dx_tot, vector<double> & dy_tot, int numNodes, int filterNum) {
    for (int i = 0; i < numNodes; i++) {            // filter and average dx and dy
        dx_filter[i].pop_front();
        dx_filter[i].push_back(dx_tot[i]);
        dy_filter[i].pop_front();
        dy_filter[i].push_back(dy_tot[i]);
    }
    for (int i = 0; i < numNodes; i++) {
        for (int j = 0; j < filterNum; j++) {
            dx_avg[i] += (dx_filter[i])[j];
            dy_avg[i] += (dy_filter[i])[j];
        }
        dx_avg[i] = dx_avg[i]/(double)filterNum;
        dy_avg[i] = dy_avg[i]/(double)filterNum;
    }
}

void computeAttractForces(SimpleGraph & graph, vector<double> & dx_tot, vector<double> & dy_tot, int edgeNum) {
    int node0 = graph.edges[edgeNum].start;
    int node1 = graph.edges[edgeNum].end;
    Node n0 = graph.nodes[node0];
    Node n1 = graph.nodes[node1];
    
    double f_attract = K_ATTRACT*((n1.x-n0.x)*(n1.x-n0.x) + (n1.y-n0.y)*(n1.y-n0.y));
    double theta = atan2(n1.y - n0.y, n1.x - n0.x);
    
    dx_tot[node0] += f_attract*cos(theta);
    dy_tot[node0] += f_attract*sin(theta);
    dx_tot[node1] -= f_attract*cos(theta);
    dy_tot[node1] -= f_attract*sin(theta);
}

void computeRepelForces(SimpleGraph & graph, vector<double> & dx_tot, vector<double> & dy_tot, int node0, int node1) {
    Node n0 = graph.nodes[node0];
    Node n1 = graph.nodes[node1];
    
    double f_repel = K_REPEL/sqrt((n1.x-n0.x)*(n1.x-n0.x) + (n1.y-n0.y)*(n1.y-n0.y));
    double theta = atan2(n1.y - n0.y, n1.x - n0.x);
    
    dx_tot[node0] -= f_repel*cos(theta);
    dy_tot[node0] -= f_repel*sin(theta);
    dx_tot[node1] += f_repel*cos(theta);
    dy_tot[node1] += f_repel*sin(theta);
}

/* Function: loadFileToGraph
 * Usage: loadFileToGraph(graph, infile)
 * -----
 * Initializes the nodes and edges in a formatted file
 */
void loadFileToGraph(SimpleGraph & graph, ifstream & file) {
    int numNodes;
    file >> numNodes;                       // read number of nodes
    initNodePositions(graph, numNodes);
    
    while (!file.fail()) {                  // read all edges
        Edge e;
        file >> e.start;
        file >> e.end;
        graph.edges.push_back(e);
    }
}

/* Function: initNodePositions
 * Usage: initNodePositions(graph, numNodes)
 * -----
 * Initializes the number of nodes specified in a circle
 */
void initNodePositions(SimpleGraph & graph, int numNodes) {
    for (int k = 0; k < numNodes; k++) {
        Node n;
        n.x = cos(2*PI*k/numNodes);
        n.y = sin(2*PI*k/numNodes);
        graph.nodes.push_back(n);
    }
}

/* Function: promptUserForLoop
 * Usage: bool loop = promptUserForLoop("Would you like to display another graph?")
 * -----
 * Returns true if the user types "yes" or "y", and false otherwise.
 */
bool promptUserForLoop(string prompt) {
    string response;
    string line;
    while (true) {
        cout << prompt;
        getline(cin, line);
        istringstream stream(line);
        stream >> response >> ws;
        if (!stream.fail() && stream.eof()) break;
        cout << "Please enter yes/no." << endl;
        if (prompt == "") prompt = "Would you like to loop?";
    }
    return (toUpperCase(response) == "YES" || toUpperCase(response) == "Y");
}

string toUpperCase(string str) {
    string result;
    for (int i = 0; i < str.size(); i++) {
        result += toupper(str[i]);
    }
    return result;
}

int promptUserForFilter(string prompt) {
    int requestedFilter;
    string line;
    while (true) {
        cout << prompt;
        getline(cin, line);
        istringstream stream(line);
        stream >> requestedFilter >> ws;
        if (!stream.fail() && stream.eof()) break;
        cout << "Illegal numeric format. Try again." << endl;
        if (prompt == "") prompt = "Enter a filter value: ";
    }
    return requestedFilter;
}

double promptUserForTime(string prompt) {
    double requestedTime;
    string line;
    while (true) {
        cout << prompt;
        getline(cin, line);
        istringstream stream(line);
        stream >> requestedTime >> ws;
        if (!stream.fail() && stream.eof()) break;
        cout << "Illegal numeric format. Try again." << endl;
        if (prompt == "") prompt = "Enter a number for time: ";
    }
    return requestedTime;
}

string promptUserForFile(string prompt) {
    string fileName;
    string line;
    while (true) {
        cout << prompt;
        getline(cin, line);
        istringstream stream(line);
        stream >> fileName >> ws;
        if (!stream.fail() && stream.eof()) break;
        cout << "Illegal file name. Try again." << endl;
        if (prompt == "") prompt = "Enter file name: ";
    }
    return fileName;
}

void Welcome() {
	cout << "Welcome to CS106L GraphViz!" << endl;
	cout << "This program uses a force-directed graph layout algorithm" << endl;
	cout << "to render sleek, snazzy pictures of various graphs." << endl;
	cout << endl;
}