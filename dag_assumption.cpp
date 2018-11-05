#include <bits/stdc++.h>
#include "carpooling.h"
using namespace std;

char * timeSourceName, * locationInputName, * edgeInputName;
char * dagOutputName, * dagDetailOutputName;
string location;

/* nodes global structure */
vector< pair<double, double> > nodesToLatLon;
vector<long long int> nodeID;
unordered_map<long long int, int> idToNode;

/* Edges global structure */
vector< vector<long long int> > edges;
vector< vector<double> > edgeWeight;
vector< vector<long long int> > edgesReverse;
vector< vector<double> > edgeWeightReverse;
vector< vector<double> > edgeTime;

vector< vector<int> > histPaths;

void getTrajectoryInput() {
	ifstream file;
	file.open( timeSourceName );
	string s;
	while( getline(file, s) ) {
		stringstream ss( s );
		int cnt = 0;
		ss >> cnt;
		vector<int> path(cnt);
		for (int i=0; i<cnt; i++) {
			ss >> path[i];
		}
		histPaths.push_back(path);
	}

	file.close();
	cout<<histPaths.size()<<endl;
}

void takeGraphInput() {
	ifstream Location;
	Location.open( locationInputName );
	string s;
	int index = 0, numNodes = 0, numEdge = 0;
	
	stringstream ss;
	char ch;
	if (location.compare("BJ") == 0) {
		getline(Location, s);
		ss.str( std::string() );
		ss.clear();
		ss<<s;
		ss>>numNodes>>ch>>numEdge;
	}

	while( getline(Location, s) ) {
		ss.str( std::string() );
		ss.clear();
		ss<<s;
		long long int id; double lon; double lat; 
		ss>>id>>ch>>lat>>ch>>lon; 
		nodesToLatLon.push_back( make_pair( lat, lon) );
		nodeID.push_back( id );
		idToNode[ id ] = index;
		// printf("%lld = %d\n",id,index);
		index++;
		if(location.compare("BJ") == 0 && nodesToLatLon.size() == numNodes)
			break;
	}

	if (location.compare("SF") == 0 || location.compare("NY") == 0 || location.compare("SG") == 0) {
		Location.close();
		Location.open( edgeInputName );
	}

	// Get edges
	int count = 0;
	edges.resize(nodeID.size());
	edgeWeight.resize(nodeID.size());
	edgeTime.resize(nodeID.size());
	edgesReverse.resize(nodeID.size());
	edgeWeightReverse.resize(nodeID.size());
	while( getline(Location, s) ) {
		ss.str( std::string() );
		ss.clear();
		ss<<s;
		long long int numRandom; long long int node1; long long int node2; double weight; char ch; int oneWay = 1; double timeNeeded;
		if (location.compare("NY") == 0 || location.compare("BJ") == 0) {
			ss>>node1>>ch>>node2>>ch>>weight>>ch>>oneWay>>ch>>timeNeeded;
		}
		else {
			ss>>node1>>ch>>node2>>ch>>weight>>ch>>oneWay;
		}
		node1 = idToNode[node1];
		node2 = idToNode[node2];
		if (location.compare("NY") == 0 || location.compare("SG") == 0 || location.compare("SF") == 0) {
			weight /= 1000;
		}
		
		// printf("edge %lld %lld\n",node1,node2);
		edges[ node1 ].push_back( node2 );
		edgeWeight[ node1 ].push_back( weight );
		// edgeTime[ node1 ].push_back( timeNeeded );
		edgeTime[ node1 ].push_back( weight / SPEED );
		edgesReverse[ node2 ].push_back( node1 );
		edgeWeightReverse[ node2 ].push_back( weight );
		count = oneWay ? count +1: count+2;
		if( !oneWay ) {
			long long int temp = node1; node1 = node2; node2 = temp;
			edges[ node1 ].push_back( node2 );
			edgeWeight[ node1 ].push_back( weight );
			// edgeTime[ node1 ].push_back( timeNeeded );
			edgeTime[ node1 ].push_back( weight / SPEED );

			edgesReverse[ node2 ].push_back( node1 );
			edgeWeightReverse[ node2 ].push_back( weight );
		}

	}
	cout<<count+1<<endl;
	Location.close();
	return ;
}

int main(int argc, char const *argv	[])
{
	if (argc < 2) {
		printf("Usage: ./a.out [location={BJ|SF|NY|SG}]\n");
		return 0;
	}

	location = argv[1];

	timeSourceName = (char*)malloc(50);
	locationInputName = (char*)malloc(50);
	edgeInputName = (char*)malloc(50);
	dagOutputName = (char*)malloc(50);
	dagDetailOutputName = (char*)malloc(50);
	if (location.compare("SF") == 0) {
		sprintf(timeSourceName, "sf_output_dag");
		sprintf(locationInputName, "sf_location");
		sprintf(edgeInputName, "sf_edge");
		sprintf(dagOutputName, "sf_dag_violate");
		sprintf(dagDetailOutputName, "sf_dag_violate_detail");
	}
	else {
		printf("Location not recognized.");
		return 1;
	}
	
	takeGraphInput();
	printf("Finish taking graph input\n"); fflush(stdout);
	getTrajectoryInput();
	printf("Finish get trajectory input\n"); fflush(stdout);

	FILE *dagOutput = fopen(dagOutputName, "w");
	FILE *dagDetailOutput = fopen(dagDetailOutputName, "w");

	for (int i=0; i<histPaths.size(); i++) {
		if (i%100 == 0) { printf("At path %d\n",i); fflush(stdout); }
		vector< double > distanceToDest( nodeID.size() );
		int source = histPaths[i][0];
		int dest = histPaths[i][histPaths[i].size()-1];
		dijkstra_lengths(nodeID.size(), dest, source, distanceToDest, edgesReverse, edgeWeightReverse);
		double now = distanceToDest[source];
		if (now == MAX_DISTANCE) continue;
		int violate = 0, prev = source;
		// printf("  %.4f\n",now);
		for (int j=1; j<histPaths[i].size(); j++) {
			int at = histPaths[i][j];
			// printf("  %.4f\n",distanceToDest[at]);
			if (distanceToDest[at] == MAX_DISTANCE) continue;
			if (distanceToDest[at] > distanceToDest[prev]) {
				fprintf(dagDetailOutput, "%.4f ",distanceToDest[at]-distanceToDest[prev]);
				violate++;
			}
			else {
				fprintf(dagDetailOutput, "%.4f ",0.0);
			}
			prev = histPaths[i][j];
		}
		fprintf(dagDetailOutput, "\n");
		fprintf(dagOutput,"%d/%d = %.4f\n",violate,(int)histPaths[i].size()-1,violate*1.0/((int)histPaths[i].size()-1));
	}
	
	return 0;
}



