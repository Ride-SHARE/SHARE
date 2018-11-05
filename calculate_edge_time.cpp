#include <bits/stdc++.h>
#include <ctime>

using namespace std;

long long n;
int queryCnt;
char * txtName, * degName, * outName, *inName;
char * timeSourceName, * locationInputName, * edgeInputName;
bool fc = 0;
string location;
int locationFactor = 1;

map< pair<int, int>, double > timeOptimize;
map< int, map<int, double> > distanceFrequentNodes;
vector< bool > frequentPickup;
vector< bool > frequentDrop;
vector< double > distanceFromSource;
vector< double > distanceFromDestination;
vector< double > distanceToDestination;

int DEBUG = 1;
bool PRINT_PATH = false;
long long int beta = 0;
long long int maxDistance = 100000;
/* nodes global structure */
vector< pair<double, double> > nodes;
vector<long long int> nodeID;
map<long long int, long long int> idToNode;
map< long long int, pair<double, double> > nodeToLatLon;

/* Edges global structure */
vector< vector<long long int> > edges;
// Dist , number of times visited
vector< vector<pair<double,int> > > edgeTime;
vector< vector<double> > edgeWeight;

vector< map<long long int, vector<double> > > sourceDestinationTime;


double distFromSourceToDestination;


void getTimeSourceDestination_Other() {
	ifstream file;
	file.open( timeSourceName );
	string s;
	int trips = 0;
	sourceDestinationTime.resize(nodeID.size());
	while( getline(file, s) ) {
		stringstream ss( s );
		long long int source, dest;
		double timeSlot;
		ss>>source>>timeSlot;
		source = idToNode[source];
		while( ss>>dest ) {
			trips++;
			dest = idToNode[dest];
			sourceDestinationTime[ source ][ dest ].push_back( timeSlot );
		}
	}
	file.close();
	cout<<trips<<endl;
}

void getTimeSourceDestination_NY() {
	ifstream file;
	file.open( timeSourceName );
	string s;
	sourceDestinationTime.resize(nodeID.size());
	int trips = 0;
	while( getline(file, s) ) {
		stringstream ss( s );
		long long int source, dest;
		double timeStart, timeEnd;
		ss>>source>>dest>>timeStart>>timeEnd;
		if (source < nodeID.size() && dest < nodeID.size()) {
			trips++;
			double minute = timeEnd - timeStart;
			while (minute < 0) {
				minute += 1440;
			}
			sourceDestinationTime[ source ][ dest ].push_back( minute );
		}

	}

	file.close();
	cout<<trips<<endl;
}

void getTimeSourceDestination() {
	if (location.compare("NY") == 0) {
		getTimeSourceDestination_NY();
	}
	else {
		getTimeSourceDestination_Other();
	}
}

void takeGraphInput() {
	ifstream Location;
	Location.open( locationInputName );
	string s;
	int index = 0, numNodes = 0, numEdge = 0;

	getline(Location, s);
	stringstream ss(s);
	char ch;
	if (location.compare("BJ") == 0) {
		ss>>numNodes>>ch>>numEdge;
	}

	while( getline(Location, s) ) {
		ss.str( std::string() );
		ss.clear();
		ss<<s;
		long long int id; double lon; double lat;
		ss>>id>>ch>>lat>>ch>>lon;
		nodes.push_back( make_pair( lon, lat) );
		nodeID.push_back( id );
		idToNode[ id ] = index;
		//printf("%lld = %d\n",id,index);
		nodeToLatLon[ id ] = make_pair( lat, lon); // not used?
		index++;
		if(location.compare("BJ") == 0 && nodes.size() == numNodes)
			break;
	}

	if (location.compare("SF") == 0 || location.compare("NY") == 0) {
		Location.close();
		Location.open( edgeInputName );
	}

	// Get edges
	int count = 0;
	edges.resize(nodeID.size());
	edgeWeight.resize(nodeID.size());
	edgeTime.resize(nodeID.size());
	while( getline(Location, s) ) {
		ss.str( std::string() );
		ss.clear();
		ss<<s;
		long long int numRandom; long long int node1; long long int node2; double weight; char ch; int oneWay = 1;
		if (location.compare("BJ") == 0) {
			ss>>numRandom>>ch>>node1>>ch>>node2>>ch>>weight>>ch>>oneWay;
		}
		else {
			ss>>node1>>ch>>node2>>ch>>weight;
		}
		node1 = idToNode[node1];
		node2 = idToNode[node2];
		if (location.compare("NY") == 0) {
			weight /= 1000;
		}

		edges[ node1 ].push_back( node2 );
		edgeWeight[ node1 ].push_back( weight );
		edgeTime[ node1 ].push_back(make_pair(0,0));
		count = oneWay ? count +1: count+2;
		if( !oneWay ) {
			long long int temp = node1; node1 = node2; node2 = temp;
			edges[ node1 ].push_back( node2 );
			edgeWeight[ node1 ].push_back( weight );
			edgeTime[ node1 ].push_back(make_pair(0,0));
		}

	}
	cout<<count+1<<endl;
	Location.close();
	return ;
}

void printEdgeOutput() {
	ofstream Location;
	Location.open( "ny_edge_time" );

	double totalmin = 0;
	int totaledge = 0;
	for (int i=0; i<edges.size(); i++) {
		for (int j=0; j<edges[i].size(); j++) {
			if (edgeTime[i][j].second > 0) {
				totalmin += (edgeTime[i][j].first / edgeTime[i][j].second) / edgeWeight[i][j];
				totaledge += 1;
			}
		}
	}

	char ch = ',';
	for (int i=0; i<edges.size(); i++) {
		double localmin = 0;
		int localedge = 0;
		for (int j=0; j<edges[i].size(); j++) {
			double et;
			if (edgeTime[i][j].second > 0) {
				et = edgeTime[i][j].first / edgeTime[i][j].second;
				//Location << "Self("<<edgeTime[i][j].first<<" "<<edgeTime[i][j].second<<") ";
				localmin += (edgeTime[i][j].first / edgeTime[i][j].second) / edgeWeight[i][j];
				localedge += 1;
			}
			else {
				continue;
			}
			Location<<nodeID[i]<<ch<<nodeID[edges[i][j]]<<ch<<edgeWeight[i][j]*1000<<ch<<'1'<<ch<<et << endl;
		}
		for (int j=0; j<edges[i].size(); j++) {
			if (edgeTime[i][j].second > 0)
				continue;
			double et;
			if (localedge > 0) {
				// Average of all edges going out from this node
				//Location << "Local ";
                et = edgeWeight[i][j] * (localmin / localedge);
			}
			else {
				// Average of global speed
				//Location << "Global ";
				et = edgeWeight[i][j] * (totalmin / totaledge);
			}
			Location<<nodeID[i]<<ch<<nodeID[edges[i][j]]<<ch<<edgeWeight[i][j]*1000<<ch<<'1'<<ch<<et << endl;
		}
	}
	Location.close();
	return ;
}

void dijkstra_lengths(long long int S, vector< double > &distanceFromSource,
	vector<long long int> &prevNode, vector< vector<long long int> > &edges, vector< vector<double> > &edgeWeight) {

	for(int i=0; i< nodeID.size(); i++)
	{ 	distanceFromSource[ i ] = 100000;
		prevNode[ i ] = -1;
	}

	distanceFromSource[ S ] = 0;
	priority_queue< pair<float, long long int> > dj;
	dj.push( make_pair(0, S) );

	pair<float, long long int> x;
	long long int u, v;
	float alt;

	while( dj.size() ) {
		x = dj.top();
		dj.pop();
		u = x.second;

		for(int i=0; i < edges[ u ].size(); i++) {
			v = edges[ u ][ i ];
			alt = distanceFromSource[ u ] + edgeWeight[ u ][ i ];
			if( alt < distanceFromSource[ v ] )
			{ 	distanceFromSource[ v ] = alt;
				dj.push( make_pair( -alt, v) );
				prevNode[ v ] = u;
			}
		}
	}
}

vector<long long int> recoverPath(long long int S, long long int D, vector<long long int> &prevNode) {
	vector<long long int> path;
	long long int node = D;
	while( true ) {
		path.push_back( node );
		if( ( node == S ) || ( node == -1) )
		  	break;
		node = prevNode[ node ];
	}

	reverse(path.begin(), path.end());

	return path;
}

void printPath(vector<long long int> &path, vector< long long int > &weights, vector< double > &distanceFromSource,	vector< double > &distanceToDestination) {
	double cumDist = 0;
	for (int i=0; i<path.size(); i++) {
		long long int v = path[i];
		if (i > 0) {
			long long int u = path[i-1];
			for (int j=0; j<edgeWeight[u].size(); j++) {
				if (edges[u][j] == v) {
					cumDist += edgeWeight[u][j];
					break;
				}
			}
		}
		printf("  %7lld:\t%lld\t(Dist: %.4f) [From source: %.4f\tTo dest: %.4f]\n",nodeID[v], weights[v], cumDist, distanceFromSource[v], distanceToDestination[v]);
	}
	fflush(stdout);
}

void addPathDist(vector<long long int> &path, double totalDist, double minuteUsed, int numTrips) {
	double cumDist = 0;
	for (int i=0; i<path.size(); i++) {
		long long int v = path[i];
		if (i > 0) {
			long long int u = path[i-1];
			for (int j=0; j<edgeWeight[u].size(); j++) {
				if (edges[u][j] == v) {
					double oldTime = edgeTime[u][j].first;
					int oldTrip = edgeTime[u][j].second;
					edgeTime[u][j] = make_pair(oldTime + minuteUsed * edgeWeight[u][j] / totalDist, oldTrip + numTrips);
					//printf("  %lld->%lld = %.4f %d\n",nodeID[u],nodeID[v],edgeTime[u][j].first,edgeTime[u][j].second);
					break;
				}
			}
		}
	}
}

void calculateEdgeTime(long long int source) {
	vector<double> distFromSource(nodeID.size());
	vector<long long int> prevNode(nodeID.size());
	//printf("Calculating\n");
	dijkstra_lengths(source, distFromSource, prevNode, edges, edgeWeight);

	for (map<long long int, vector<double> >::iterator iter = sourceDestinationTime[source].begin(); iter != sourceDestinationTime[source].end(); iter++) {
		vector<double> timeList = iter->second;
		double minuteUsed = 0;
		int numTrips = timeList.size();
		for (int i=0; i<numTrips; i++) {
			minuteUsed += timeList[i];
		}

		long long int destination = iter->first;
		vector<long long int> path = recoverPath(source, destination, prevNode);
		double totalDist = distFromSource[destination];
		if(totalDist > 99999 ) {
			//printf("Cannot reach destination?\n");
			continue;
		}

		addPathDist(path, totalDist, minuteUsed, numTrips);
		//printf("Done = %.4f (%.4f, %d)\n",totalDist, minuteUsed, numTrips);
	}

}

/* edges global structure */
int main(int argc, char const *argv	[])
{
	if (argc < 2) {
		printf("Usage: ./a.out [location={BJ|SF|NY}]\n");
		return 0;
	}
	location = argv[1];
	txtName = (char*)malloc(50);
	timeSourceName = (char*)malloc(50);
	locationInputName = (char*)malloc(50);
	edgeInputName = (char*)malloc(50);
	if (location.compare("BJ") == 0) {
		sprintf(txtName, "beijingIndex");
		sprintf(timeSourceName, "bj_output");
		sprintf(locationInputName, "bj_graph");
	}
	else if (location.compare("NY") == 0) {
		sprintf(txtName, "nyIndex");
		sprintf(timeSourceName, "ny_output_timestamp");
		sprintf(locationInputName, "ny_location");
		sprintf(edgeInputName, "ny_edge");
		locationFactor = 2;
	}
	else if (location.compare("SF") == 0) {
		sprintf(txtName, "sfIndex");
		sprintf(timeSourceName, "sf_output");
		sprintf(locationInputName, "sf_location");
		sprintf(edgeInputName, "sf_edge");
	}
	else {
		printf("Location not recognized.");
		return 1;
	}

	takeGraphInput();
	printf("Finish taking graph input\n"); fflush(stdout);
	getTimeSourceDestination();
	printf("Finish get time source destination\n"); fflush(stdout);

	for( int i=0; i < sourceDestinationTime.size(); i += 1) {
		if (i%10 == 0) printf("i = %d\n",i);
		calculateEdgeTime(i);
	}

	printEdgeOutput();

	printf("END OF SIMULATION\n");

	return 0;
}

