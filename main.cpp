#include <bits/stdc++.h>
#include <ctime>
#include "carpooling.h"

using namespace std;

char * timeSourceName, * locationInputName, * edgeInputName;
char * txtName;
string location;

/* Parameters */
int DEBUG = 1;
bool PRINT_PATH = false;
bool FREQUENCY_OPT = false;
bool SUBSET = true;
int TOTAL_TRIP = 50;
double alpha  = 1.3;
double maxDepth = 0.2;
long long int maxDistance = 100000;
/* nodes global structure */
vector< pair<double, double> > nodesToLatLon;
vector<long long int> nodeID;
map<long long int, long long int> idToNode;

/* Edges global structure */
vector< vector<long long int> > edges;
vector< vector<double> > edgeWeight;
vector< vector<long long int> > edgesReverse;
vector< vector<double> > edgeWeightReverse;

vector< map<long long int, vector<long long int> > > sourceTimeDestination; 
vector< map<long long int, vector<long long int> > > simulation; 

int dijkstraScore[1000];
double dijkstraTime[1000];
double dijkstraDist[1000];

int maxScorePerDistanceScore[1000];
double maxScorePerDistanceTime[1000] ;
double maxScorePerDistanceDist[1000];

int dagScore[1000];
double dagTimeTaken[1000];
double dagDist[1000];

int dagExScore[1000];
double dagExTime[1000];
double dagExDist[1000];

long long int tripNumber = 0;
clock_t startTimeTrip;
clock_t endTimeTrip;

void getTimeSourceDestination_Other() {
	ifstream file;
	file.open( timeSourceName );
	string s;
	int trips = 0;
	sourceTimeDestination.resize(nodeID.size());
	simulation.resize(nodeID.size());
	while( getline(file, s) ) {
		stringstream ss( s );
		long long int source, timeSlot, dest;
		ss>>source>>timeSlot;
		source = idToNode[source];
		while( ss>>dest ) {
			trips++;
			dest = idToNode[dest];
			int prob = rand()%100;
			if( prob < 80 )
				sourceTimeDestination[ source ][ timeSlot ].push_back( dest );
			else 
				simulation[ source ][ timeSlot ].push_back( dest );
		}
	}
	file.close();
	cout<<trips<<endl;
}

void getTimeSourceDestination_NY() {
	ifstream file;
	file.open( timeSourceName );
	string s;
	sourceTimeDestination.resize(nodeID.size());
	simulation.resize(nodeID.size());
	int trips = 0;
	while( getline(file, s) ) {
		stringstream ss( s );
		long long int source, timeSlot, dest;
		ss>>source>>timeSlot>>dest;
		if (source < nodeID.size() && dest < nodeID.size()) {
			trips++;

			int prob = rand()%100;
			if( prob < 80 )
				sourceTimeDestination[ source ][ timeSlot ].push_back( dest );
			else 
				simulation[ source ][ timeSlot ].push_back( dest );			
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
		nodesToLatLon.push_back( make_pair( lat, lon) );
		nodeID.push_back( id );
		idToNode[ id ] = index;
		//printf("%lld = %d\n",id,index);
		index++;
		if(location.compare("BJ") == 0 && nodesToLatLon.size() == numNodes)
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
	edgesReverse.resize(nodeID.size());
	edgeWeightReverse.resize(nodeID.size());
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
		edgesReverse[ node2 ].push_back( node1 );
		edgeWeightReverse[ node2 ].push_back( weight );
		count = oneWay ? count +1: count+2;
		if( !oneWay ) {
			long long int temp = node1; node1 = node2; node2 = temp;
			edges[ node1 ].push_back( node2 );
			edgeWeight[ node1 ].push_back( weight );

			edgesReverse[ node2 ].push_back( node1 );
			edgeWeightReverse[ node2 ].push_back( weight );
		}

	}
	cout<<count+1<<endl;
	Location.close();
	return ;
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

int getPathScore(vector<long long int> &path, vector< long long int > &weights) {
	int score = 0;
	set<long long int> pathNode;
	for(int i=1;i<path.size(); i++) {
		if (pathNode.find(path[i]) != pathNode.end()) {
			// This path has cycle
			printf("ERROR: Cycle exists at node %lld\n",nodeID[path[i]]);
			return -1;
		}
		score += weights[ path[i] ];
		pathNode.insert(path[i]);
	}
	return score;
}

double getPathDist(vector<long long int> &path) {
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
	}
	return cumDist;
}


bool potentialScore(long long int source, long long int destination, long long int timeSlot)
{
	cout<<"Trip start"<<endl;
	printf("Source: %lld(%lld)\tDestination: %lld(%lld)\tTime slot: %lld\n",source,nodeID[source],destination,nodeID[destination],timeSlot);

	vector< double > distanceFromSource(nodeID.size());
	vector< double > distanceFromDestination(nodeID.size());
	vector< double > distanceToDestination(nodeID.size());

	startTimeTrip = clock();
	vector<long long int> dijkstraPath = dijkstra_lengths(nodeID.size(), source, destination, distanceFromSource, edges, edgeWeight);
	dijkstraTime[ tripNumber ] = float( clock()- startTimeTrip )/ CLOCKS_PER_SEC ;

	if(distanceFromSource[destination] == maxDistance || !distanceFromSource[ destination ] ) {
		printf("Cannot reach destination?\n");
		return false;
	}

	dijkstra_lengths(nodeID.size(), destination, source, distanceToDestination, edgesReverse, edgeWeightReverse);
	double timeDij = ( clock() - startTimeTrip ) / (float) CLOCKS_PER_SEC;

	dijkstra_lengths(nodeID.size(), destination, source, distanceFromDestination, edges, edgeWeight);

	if( DEBUG ) {
		printf("  [Dijkstra length] %.4f\n",distanceFromSource[ destination ]);
		printf("  [Query length] %.4f\n",query(source, destination));
	}

	startTimeTrip = clock();

	vector<long long int> expectedTrips( nodeID.size() ); 
	// From s to d, via v, where v is within 1+alpha
	for( int i =0; i< nodeID.size(); i++) {	
		if( ( distanceFromSource[ i ] + distanceToDestination[ i ] ) 
				< alpha	*distanceFromSource[ destination ] ) {
			//long long int passengersNumber;
			expectedTrips[ i ] = get_expected_trips(source, destination, i, timeSlot, alpha, sourceTimeDestination, distanceFromSource, distanceToDestination, distanceFromDestination);
		}
	}

	vector<long long int> simulationDataset( nodeID.size() ); 
	// From s to d, via v, where v is within 1+alpha
	for( int i =0; i<nodeID.size(); i++) {	
		if( ( distanceFromSource[ i ] + distanceToDestination[ i ] ) 
				< alpha	*distanceFromSource[ destination ] ) {
			//long long int passengersNumber;
			simulationDataset[ i ] = get_expected_trips(source, destination, i, timeSlot, alpha, simulation, distanceFromSource, distanceToDestination, distanceFromDestination);
		}
	}
	cout<<"TIME: "<<( clock() - startTimeTrip ) / (float) CLOCKS_PER_SEC<<" TIME COMPLETE\n";
	/* Dijkstra path */
	int dijkstraPathScore = getPathScore(dijkstraPath, expectedTrips);
	dijkstraScore[ tripNumber ] = dijkstraPathScore+1;
	dijkstraDist[ tripNumber ] = getPathDist(dijkstraPath);
	printf("  [Dijkstra score] %d (Time: %.2f)\n",dijkstraPathScore  + 1,dijkstraTime[ tripNumber ]);
	if (PRINT_PATH) {
		printPath(dijkstraPath, expectedTrips, distanceFromSource, distanceToDestination);
	}

	/* Max Score / Distance path */
	startTimeTrip = clock();
	
	vector<long long int> maxScorePerDistancePath = greedyScoreDivDist( nodeID.size(), source, destination, alpha, expectedTrips, 
		distanceFromSource, edges, edgeWeight);
	int maxScorePerDistancePathScore = getPathScore(maxScorePerDistancePath, expectedTrips);
	maxScorePerDistanceTime[ tripNumber ] = ( clock() - startTimeTrip ) / (float) CLOCKS_PER_SEC + timeDij;
	maxScorePerDistanceScore[ tripNumber ] = maxScorePerDistancePathScore+1;
	maxScorePerDistanceDist[ tripNumber ] = getPathDist(maxScorePerDistancePath);
	printf("  [Max Score/Distance path score] %d (Time: %.2f)\n",maxScorePerDistancePathScore + 1, 
						maxScorePerDistanceTime[ tripNumber ] );
	if (PRINT_PATH) {
		printPath(maxScorePerDistancePath, expectedTrips, distanceFromSource, distanceToDestination);
	}

	/* calculating DAG parameters */
	startTimeTrip = clock();
	vector<long long int> dagPath = findDAGPath(nodeID.size(), source, destination, timeSlot, alpha, sourceTimeDestination,
		distanceFromSource,	distanceToDestination, distanceFromDestination, edges, edgeWeight);
	int scoreDag = getPathScore(dagPath, expectedTrips);

	dagTimeTaken[ tripNumber ] = float( clock()- startTimeTrip )/ CLOCKS_PER_SEC + timeDij;
	printf("  [DAG score] %d (Time: %.2f)\n",scoreDag+1,dagTimeTaken[ tripNumber ]);
	if (PRINT_PATH) {
		printPath(dagPath, expectedTrips, distanceFromSource, distanceToDestination);
	}
	dagScore[ tripNumber ] = scoreDag + 1 ;
	dagDist[ tripNumber ] = getPathDist(dagPath);
	
	/* calculating back DAG	parameters */
	startTimeTrip = clock();

	vector<long long int> dagExPath = findDAGExtendedPath( nodeID.size(), source, destination, timeSlot, alpha, maxDepth, sourceTimeDestination,
										 distanceFromSource, distanceToDestination, distanceFromDestination);
	dagExTime[ tripNumber ] = ( clock() - startTimeTrip ) / (float) CLOCKS_PER_SEC + timeDij;
	int scoreExDag = getPathScore(dagExPath, expectedTrips);
	dagExScore[ tripNumber ] = scoreExDag + 1 ;
	dagExDist[ tripNumber ] = getPathDist(dagExPath);
	printf("  [Extended DAG score] %d (Time: %.2f)\n",scoreExDag+1,dagExTime[ tripNumber ]);
	if (PRINT_PATH) {
		printPath(dagExPath, expectedTrips, distanceFromSource, distanceToDestination);
	}

	// remove cached data
	timeOptimize.clear();

	return true;	
}


/* edges global structure */
int main(int argc, char const *argv	[])
{
	if (argc < 4) {
		printf("Usage: ./a.out [location={BJ|SF|NY}] [maxDepth=0.2] [alpha=1.3]\n");
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
		sprintf(timeSourceName, "ny_output");
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
	
	string s = argv[2];
	stringstream ss(s);	
	stringstream ssDebug(argv[3]);
	ss>>maxDepth;
	ssDebug>>alpha;
	printf("Max Depth = %.2f\nAlpha = %.2f\n",maxDepth,alpha); fflush(stdout);

	takeGraphInput();
	printf("Finish taking graph input\n"); fflush(stdout);
	getTimeSourceDestination();
	printf("Finish get time source destination\n"); fflush(stdout);

	if (SUBSET) {
		// Map old index to new index
		map<int, int> marked;
		vector<long long int> nodeIDTrunc;
		for( int i =0 ; i < nodeID.size() ; i++) {
			if( nodesToLatLon[ i ].first < 40.6097  ) {
				marked[i] = nodeIDTrunc.size();
				nodeIDTrunc.push_back(nodeID[i]);
			}
		}
		nodeID = nodeIDTrunc;
		printf("Number of marked node: %d\n",(int)nodeID.size()); fflush(stdout);
		// nodeID, edges, edgeWeight, edgesReverse, edgeWeightReverse, sourceTimeDestination, simulation
		subset_edge(marked, edges, edgeWeight);
		subset_edge(marked, edgesReverse, edgeWeightReverse);
		subset_dataset(marked, sourceTimeDestination);
		subset_dataset(marked, simulation);
	}


	queryInit(txtName);
	intializeFrequent(nodeID.size(), FREQUENCY_OPT, sourceTimeDestination, edges, edgeWeight);
	
	assignExtendEdge(nodeID.size(), maxDepth, edges, edgeWeight);

	//srand (time(NULL));
	srand(0);
	
	tripNumber = 0;
	for( int i=0; tripNumber < TOTAL_TRIP; i += 1) {
		long long int source, timeS, destination;
			
		source = rand() % nodeID.size() ;
		timeS = rand()%96;
		
		if( source >= sourceTimeDestination.size() )
			continue;
		else if ( sourceTimeDestination[ source ].find( timeS ) == sourceTimeDestination[ source ].end() )
			continue;
		else if( !sourceTimeDestination[ source ][ timeS ].size() )
			continue;

		destination = sourceTimeDestination[ source ][ timeS ][ rand() % sourceTimeDestination[ source ][ timeS ].size()  ];
		bool reached = potentialScore(source, destination, timeS);
		if (!reached) {
			continue;
		}

		printf("Trip\t%lld\n",tripNumber);
		printf("Score\t%d\t%d\t%d\t%d\n",dijkstraScore[tripNumber],maxScorePerDistanceScore[tripNumber],dagScore[tripNumber],dagExScore[tripNumber]);
		printf("Time\t%.4f\t%.4f\t%.4f\t%.4f\n",dijkstraTime[tripNumber],maxScorePerDistanceTime[tripNumber],dagTimeTaken[tripNumber],dagExTime[tripNumber]);
		printf("Dist\t%.4f\t%.4f\t%.4f\t%.4f\n",dijkstraDist[tripNumber],maxScorePerDistanceDist[tripNumber],dagDist[tripNumber],dagExDist[tripNumber]);
		printf("\n");
		tripNumber++;

	}

	printf("END OF SIMULATION\n");

	return 0;
}

