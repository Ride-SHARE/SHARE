#include <bits/stdc++.h>
#include <ctime>
#include "carpooling.h"

using namespace std;

char * timeSourceName, * locationInputName, * edgeInputName;
char * txtName;
string location;
string algoUsed;

/* nodes global structure */
vector< pair<double, double> > nodesToLatLon;
vector<long long int> nodeID;
map<long long int, long long int> idToNode;

/* Edges global structure */
vector< vector<double> > edgeTime;
vector< vector<long long int> > edgesReverse;
vector< vector<double> > edgeWeightReverse;

vector< map<long long int, vector<long long int> > > sourceTimeDestination; 
vector< map<long long int, vector<long long int> > > simulation;

vector< pair<int, pair<long long int, long long int> > > sample;

int dijkstraScore[1000];
double dijkstraTime[1000];
double dijkstraDist[1000];

int dagScore[1000];
double dagTimeTaken[1000];
double dagDist[1000];

int dagExScore[1000];
double dagExTime[1000];
double dagExDist[1000];

int dagExScore2[1000];
double dagExTime2[1000];
double dagExDist2[1000];

int tripNumber = 0;
/* Simulation */
int DAY = 6;
int TOTAL_TRIP = 200;
/*
* MODE 1: Infinity cab. After passengerWaitFactor timeSlot, the passenger will be picked up immediately
* MODE 2: Limited cab (availableCab). Infinity passengerWaitFactor
*/

vector< passenger > passengerPickedList;
// At source v - a list of passenger waiting for cab
map<long long int, vector< passenger > > passengerQueue;


int getPathScoreModified(vector<long long int> &path, int timeSlot, 
	vector< map<long long int, vector<long long int> > > &dataset,
 	const vector< double > &distanceFromSource, const vector< double > &distanceToDestination, const vector< double > &distanceFromDestination ) {
	int score = 0;
	double edgeDist = 0;

	for(int i = 1;i < path.size() - 1; i++) {
		bool found = false;
		for( int j = i + 1; j < path.size()-1; j++) {
			if( path[ i ] == path[ j ] ) {
				cout<<"CYCLE DETECTED"<<endl;
				break;
			}
		}

		if( (path[ i ] < 0) || (path[ i ] >= edges.size()) )
			continue;

		for( int edgeIndex = 0; edgeIndex < edges[ path [ i ] ].size(); edgeIndex++){
			if( edges[ path[ i ] ][ edgeIndex ] == path[ i + 1 ] ) {
				found = true;
				edgeDist += edgeWeight[ path[ i ] ][ edgeIndex ];	
				break;
			}
		}
		if( !found ) {
			return -2;
		}

		score += get_expected_trips_new( path[ 0 ], path[ path.size() -1 ], path[ i ], 
		edgeDist, timeSlot + int( edgeDist/(SPEED*deltaTime) ), alpha, dataset,
		distanceFromSource, distanceToDestination, distanceFromDestination);
	}

	return score;
}

vector<long long int> getBestPath(long long int source, long long int destination, long long int timeSlot, int reachOnly = 0, double alphaChanged = alpha ) {
	cout<<"Trip start"<<endl;
	int start_time = clock();
	printf("Source: %lld(%lld)\tDestination: %lld(%lld)\tTime slot: %lld\n",source,nodeID[source],destination,nodeID[destination],timeSlot); fflush(stdout);

	vector< double > distanceFromSource( nodeID.size() );
	vector< double > distanceFromDestination( nodeID.size() );
	vector< double > distanceToDestination( nodeID.size() );

	vector<long long int> dagExPath;
	dagExPath = dijkstra_lengths(nodeID.size(), source, destination, distanceFromSource, edges, edgeWeight);

	if( (distanceFromSource[ destination ] == maxDistance )  || !distanceFromSource[ destination ] ) {
		printf("Cannot reach destination?\n");
		vector<long long int> emptyPath;
		return emptyPath;
	}
	
	if (reachOnly) {
		vector<long long int> emptyPath(1);
		return emptyPath;
	}

	dijkstra_lengths(nodeID.size(), destination, source, distanceToDestination, edgesReverse, edgeWeightReverse);
	dijkstra_lengths(nodeID.size(), destination, source, distanceFromDestination, edges, edgeWeight);
	

	if (ALGO == 1) {
		dagExPath = findDAGPath( nodeID.size(), source, destination, timeSlot/3, alphaChanged, sourceTimeDestination, distanceFromSource, distanceToDestination, distanceFromDestination, edges, edgeWeight);
	}
	else if (ALGO == 2) {
		printf("  Start DAG EX\n"); fflush(stdout);
		dagExPath = findDAGExtendedPath( nodeID.size(), source, destination, timeSlot/15, alphaChanged, maxDepth, sourceTimeDestination, distanceFromSource, distanceToDestination, distanceFromDestination);
	
	}
	else if (ALGO == 3) {
		printf(" hueristicNew\n"); fflush(stdout);
		dagExPath = hueristicNew( nodeID.size(), source, destination, timeSlot/15, alphaChanged, sourceTimeDestination, distanceFromSource, distanceToDestination, distanceFromDestination, edges, edgeWeight);
	
	}

	if( !ALGO ) {
		dijkstraTime[ tripNumber ] = float( clock( ) - start_time)/(float) CLOCKS_PER_SEC;
		dijkstraScore[ tripNumber ] = getPathScoreModified( dagExPath, timeSlot, simulation, distanceFromSource, distanceToDestination, distanceFromDestination);
		dijkstraDist[ tripNumber ] = checkAlphaPath( dagExPath, edges, edgeWeight,
										distanceFromSource[ destination ], alpha);
	} 
	else if( ALGO == 1 ) {
		dagTimeTaken[ tripNumber ] = float( clock( ) - start_time)/(float) CLOCKS_PER_SEC;
		dagScore[ tripNumber ]	= getPathScoreModified( dagExPath, timeSlot, simulation, distanceFromSource, distanceToDestination, distanceFromDestination);
		dagDist[ tripNumber ] = checkAlphaPath( dagExPath, edges, edgeWeight,
										distanceFromSource[ destination ], alpha);
	}
	else if ( ALGO == 2 ){
		dagExTime[ tripNumber ]  = float( clock( ) - start_time)/(float) CLOCKS_PER_SEC;
		dagExScore[ tripNumber ] = getPathScoreModified( dagExPath, timeSlot, simulation, distanceFromSource, distanceToDestination, distanceFromDestination);
		dagExDist[ tripNumber ] = checkAlphaPath( dagExPath, edges, edgeWeight,
										distanceFromSource[ destination ], alpha);
	}
	else if ( ALGO == 3 ){
		dagExTime2[ tripNumber ]  = float( clock( ) - start_time)/(float) CLOCKS_PER_SEC;
		dagExScore2[ tripNumber ] = getPathScoreModified( dagExPath, timeSlot, simulation, distanceFromSource, distanceToDestination, distanceFromDestination);
		dagExDist2[ tripNumber ] = checkAlphaPath( dagExPath, edges, edgeWeight,
										distanceFromSource[ destination ], alpha);
	}

	cout<<"Trip Found"<<endl; fflush(stdout);
	
	// cout<<"CHECK ALPHA VIOLATION AMOUNT: "<<checkAlphaPath( dagExPath, edges, edgeWeight,
	// distanceFromSource[ destination ], alpha)<<endl;

	dagExPath.push_back(1);
	return dagExPath;
}


void getTimeSourceDestination_Other() {
	ifstream file;
	file.open( timeSourceName );
	string s;
	int trips = 0;
	sourceTimeDestination.resize(nodeID.size());
	simulation.resize( maxEndTime / deltaTime );
	while( getline(file, s) ) {
		stringstream ss( s );
		long long int date, source, timeSlot, dest;
		ss>>date>>source>>timeSlot;
		if (date > DAY) continue;
		if (!(timeSlot >= startTime / deltaTime)) continue;
		source = idToNode[source];
		while( ss>>dest ) {
			trips++;
			dest = idToNode[dest];
			if (date < DAY)
				sourceTimeDestination[ source ][ timeSlot ].push_back( dest );
			else {
				if (timeSlot < endTime / deltaTime)
					simulation[ timeSlot ][ source ].push_back( dest );
			}
		}
	}
	file.close();
	cout<<trips<<endl;
}

void getTimeSourceDestination_NYSG() {
	ifstream file;
	file.open( timeSourceName );
	string s;
	sourceTimeDestination.resize( nodeID.size() );
	simulation.resize( nodeID.size() );
	int trips = 0;
	while( getline(file, s) ) {
		stringstream ss( s );
		long long int date, source, timeSlot, dest;
		ss>>date>>source>>timeSlot>>dest;
		if (location.compare("SG") == 0) {
			timeSlot /= deltaTime;
		}
		//printf("Hello %lld %lld %lld %lld\n",date,source,timeSlot,dest);
		if (source < nodeID.size() && dest < nodeID.size() && timeSlot >= startTime / deltaTime) {
			trips++;

			if( date < DAY )
				sourceTimeDestination[ source ][ timeSlot/15 ].push_back( dest );
			else {
				if (timeSlot < endTime / deltaTime) {
					simulation[ source ][ timeSlot/5 ].push_back( dest );		
					sample.push_back( make_pair( timeSlot/5, make_pair(source, dest) ) ); 	
				}
			}
		}
		
	}

	file.close();
	cout<<trips<<endl;
}

void getTimeSourceDestination() {
	if (location.compare("NY") == 0 || location.compare("SG") == 0) {
		getTimeSourceDestination_NYSG();
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
		if (location.compare("SG") == 0) {
			ss>>node1>>ch>>node2>>ch>>weight>>ch>>oneWay;
		}
		else {
			ss>>node1>>ch>>node2>>ch>>weight>>ch>>oneWay>>ch>>timeNeeded;
		}
		node1 = idToNode[node1];
		node2 = idToNode[node2];
		if (location.compare("NY") == 0 || location.compare("SG") == 0) {
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

vector< long long int > potentialScore(long long int source, long long int destination, long long int timeSlot)
{
	vector< long long int > path;
	/* Max Score / Distance path */
	for( int i = 0; i < 4; i++) { 
		ALGO = i;
		path = getBestPath( source, destination, timeSlot);
	}

	return path;	
}

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
		sprintf(timeSourceName, "bj_output_simulation_5");
		sprintf(locationInputName, "bj_graph_time");
		DAY = 6;
	}
	else if (location.compare("NY") == 0) {
		sprintf(txtName, "nyIndex");
		sprintf(timeSourceName, "ny_output");
		sprintf(locationInputName, "ny_location");
		sprintf(edgeInputName, "ny_edge_time");
		DAY = 24;
		// locationFactor = 2;
	}
	else if (location.compare("SF") == 0) {
		sprintf(txtName, "sfIndex");
		sprintf(timeSourceName, "sf_output");
		sprintf(locationInputName, "sf_location");
		sprintf(edgeInputName, "sf_edge");
	}
	else if (location.compare("SG") == 0) {
		sprintf(txtName, "sgIndex");
		sprintf(timeSourceName, "sg_output");
		sprintf(locationInputName, "sg_location");
		sprintf(edgeInputName, "sg_edge");
		DAY = 16;
	}
	else {
		printf("Location not recognized.");
		return 1;
	}
	

	degName = (char*)malloc(1+strlen(txtName) + 50);
	sprintf(degName, "%s.deg", txtName);
	
	inName = (char*)malloc(1+strlen(txtName) + 50);
	sprintf(inName, "%s.labelin", txtName);
	outName = (char*)malloc(1+strlen(txtName) + 50);
	sprintf(outName, "%s.labelout", txtName);

	timer tm;
	
	loadIndex();
	
	printf("load time %lf (ms)\n", tm.getTime()*1000); fflush(stdout);
	
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

	queryInit( txtName );
	intializeFrequent(nodeID.size(), FREQUENCY_OPT, sourceTimeDestination, edges, edgeWeight);

	assignExtendEdge(nodeID.size(), maxDepth, edges, edgeWeight);
	srand( 0 );
	for( int i = 0; tripNumber < TOTAL_TRIP; i += 1) {
		cout<<tripNumber<<endl;
		int index = rand() % sample.size() ;
		bool val = false;
		// if( query(sample[ index ].second.first, sample[ index ].second.second)  )
		// 	continue;
		val = potentialScore(sample[ index ].second.first, sample[ index ].second.second, sample[ index ].first).size() ? true: false;
		
		printf("Trip\t%lld\n",tripNumber);
		printf("Score\t%d\t%d\t%d\t%d\n",dijkstraScore[tripNumber], dagScore[tripNumber], dagExScore[tripNumber], dagExScore2[tripNumber]);
		printf("Time\t%.4f\t%.4f\t%.4f\t%.4f\n",dijkstraTime[tripNumber], dagTimeTaken[tripNumber], dagExTime[tripNumber], dagExTime2[tripNumber]);
		printf("Detour\t%.4f\t%.4f\t%.4f\t%.4f\n",dijkstraDist[tripNumber], dagDist[tripNumber], dagExDist[tripNumber], dagExDist2[tripNumber]);
		printf("\n");

		if( val ) {
			tripNumber++;
		}
	}

	return 0;
}

