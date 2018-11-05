#include <bits/stdc++.h>
#include "ioD.h"
#include "Hungarian.h"
using namespace std;
// #define unordered_set set
// #define unordered_map map

/* Constants */
int DEBUG = 0;
bool PRINT_PATH = false;
bool FREQUENCY_OPT = false;

const int MAX_CAB_CAPACITY = 3;
const double SPEED = 1; // km per min
const int MINUTE_PER_HISTORY_SLOT = 15;
const int MINUTE_BEFORE_ASSIGN_EMPTY_CAB = 5;
const int DELTA_TIME = 1;
const int MAX_END_TIME = 24*60;
const int MAX_DISTANCE = 1000000;
const int MINUTE_BEFORE_REJECTION = 15;
const int DISCRETE_TIME = 100;
const double DELTA_DETOUR = 0.2;

const bool TIME_VARYING = false;
const bool USE_DECAY= false;
const bool RUN_OPTIMAL= false;
const bool USE_CLUSTER= true;

// For price aware comparison only
static double BCOST;
static double FUEL;


/* Parameters */
int ASSIGN = 1;
int ROUTE = 0;
int availableCab = 800;
int cabCapacity = 2;
double alpha  = 1.3;
double maxDepth = 0.2;
int startTime  = 0;
int endTime = 24*60;

struct request {
	int id;
	int timeSlot;
	int source;
	int destination;
	double revenue;
};

struct passenger {
	int id;
	int startTimeSlot;
	int pickedTimeSlot;
	int source;
	int destination;
	int passengerNo; // First passenger / Second passenger
	bool iShared;
	double shortestPathDist;
	double actualTravelDist;
	double revenue;
};

struct event {
	int pid;
	int node;
	bool isDropoff;
	event(int apid, int anode, bool aisDropoff) : pid(apid), node(anode), isDropoff(aisDropoff) {}
};

// Query utilities

char * degName, * outName, *inName;
edgeL * deg;
edgeS * labelout, *labelin;
edgeS * labelx, * labely;

map< pair<int, int>, double > timeOptimize; 
map< int, map<int, double> > distanceFrequentNodes;
vector< bool > frequentPickup;
vector< bool > frequentDrop;

// Takes in node index
double query(int x, int y)
{	
	if (x == y) return 0;

	// if we already have this key pair value then return 
	if( frequentPickup[ x ] && frequentDrop[ y ] ) {
		return distanceFrequentNodes[ x ][ y ];
	}
	
	if( timeOptimize.find( make_pair(x, y) ) != timeOptimize.end() ) {
		return timeOptimize[ make_pair(x, y) ]; 
	}

	int xx = x, yy = y;

	x = ((deg[xx].x<<32)>>32);
	y = ((deg[yy].x<<32)>>32);
		
	if (x > y)
	{
		labelx = labelout + deg[xx].w;
		labely = labelin + deg[yy].y;
	}
	else
	{
		int xy = x; x = y; y = xy;
		labelx = labelin + deg[yy].y;
		labely = labelout + deg[xx].w;
	}

	int ans = 1000000, i = 0, j = 0;

	if (labelx[i].x != -1 && labely[j].x != -1)
	while (labelx[i].x < y)
	{
		if (labelx[i].x == labely[j].x) 
		{
			ans = ans>(labelx[i].w + labely[j].w)?(labelx[i].w + labely[j].w):ans;
			if (labelx[++i].x == -1) break;
			if (labely[++j].x == -1) break;
		}
		else if (labelx[i].x < labely[j].x)
		{
			if (labelx[++i].x == -1) break;
		}
		else if (labely[++j].x == -1) break;
	}
	
	while (labelx[i].x != -1 && labelx[i].x < y) i++;
	if (labelx[i].x == y) ans = ans>labelx[i].w?labelx[i].w:ans;

	// save the key-pair value here 
	// Note that x and y are changed during calculation, we have to use xx and yy to store the timeOptimize
	timeOptimize[ make_pair(xx, yy) ] = float(ans)/1000;
	
	return float(ans)/1000;
}

void loadIndex()
{
	long long n;
	inBufL degBuf(degName);
	inBufS inLabel(inName), outLabel(outName);
	
	n = checkB(degName)/sizeof(edgeL);

	deg = (edgeL *)malloc(sizeof(edgeL)*n);
	labelin = (edgeS*)malloc(checkB(inName));
	labelout = (edgeS*)malloc(checkB(outName));

	printf("%lld vertices\n", n);

	degBuf.start();
	for (int i = 0; i < n; i++)
		degBuf.nextEdge(deg[i]);

	inLabel.start();
	for (int i = 0; !inLabel.isEnd; i++)
		inLabel.nextEdge(labelin[i]);
	
	outLabel.start();
	for (int i = 0; !outLabel.isEnd; i++)
		outLabel.nextEdge(labelout[i]);			
}

int getPathScore(vector<int> &path, vector< long long int > &weights) {
	int score = 0;
	for(int i = 1;i < path.size() - 1; i++) {
		score += weights[ path[ i ] ];
	}
	return score;
}

double getPathDist(vector<int> &path, const vector< vector< long long int> > &edges, const vector< vector<double> > &edgeWeight) {
	double cumDist = 0;
	for (int i=0; i<path.size(); i++) {
		int v = path[i];
		if (i > 0) {
			int u = path[i-1];
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

// Check Alpha constraint for path output
double checkAlphaPath( vector< int > &path, vector< vector< long long int> > &edges, vector< vector<double> > &edgeWeight, double shortDistance) {	
	double pathLen = getPathDist(path, edges, edgeWeight);
	return (pathLen / shortDistance);
}


vector<int> dijkstra_lengths(int N, int S, int D, vector< double > &distanceFromSource,
	const vector< vector<long long int> > &edges,const vector< vector<double> > &edgeWeight);

/* new function: we store distance to and from nodes which are in 
top 10 percentile wrt number of trips from them. 
*/
void intializeFrequent(int N, bool frequentOptimize, vector< unordered_map<int, vector<request> > > &sourceTimeDestination,
	vector< vector<long long int> > &edges, vector< vector<double> > &edgeWeight) {

	frequentPickup.resize( N, false);
	frequentDrop.resize( N, false);
	if (!frequentOptimize) {
		return;
	}

	printf("INITIALIZE SHORTEST DISTANCE FOR FREQUENT NODES\n"); fflush(stdout);

	// frequent pickups
	vector<int> v;
	for( int i = 0; i< N; i++ ) {
		int count = 0;
		for( int j = 0; j < 96; j++) {
			count += sourceTimeDestination[ i ][ j ].size() ; 
		}
		v.push_back( count );
	}
	sort( v.begin(), v.end() );
	int threshold = v[ v.size()*90/100 ];
	
	for( int i = 0; i< N; i++ ) {
		int count = 0;
		for( int j = 0; j < 96; j++) {
			count += sourceTimeDestination[ i ][ j ].size() ; 
		}
		
		if(  count >= threshold ) {
			frequentPickup[ i ] = 1;
		}
		else 
			frequentPickup[ i ] = 0; 
	}

	// frequent drops offs
	vector<int> des( N, 0);
	for( int i = 0; i< N; i++ ) {
		for( int j = 0; j < 96; j++) {
			for( int k = 0; k < sourceTimeDestination[ i ][ j ].size(); k++) {
				des[ sourceTimeDestination[ i ][ j ][ k ].destination ] += 1;
			} 
		}
	}

	vector< int > desSort = des; 
	sort( desSort.begin(), desSort.end() );
	threshold = desSort[ desSort.size()*90/100 ];
	
	for( int i = 0; i< N; i++ ) {
		if(  des[ i ] >= threshold ) {
			frequentDrop[ i ] = 1;
		}
		else 
			frequentDrop[ i ] = 0; 
	}

	for( int i = 0; i< N; i++ ) {	
		if( !frequentPickup[ i ] ) 
			continue;

	//	cout<<i<<" "<<nodeID.size()<<endl;
		vector< double > distanceFromSourceL( N );
		dijkstra_lengths(N, i, -1, distanceFromSourceL, edges, edgeWeight);
	
		for( int j = 0;  j < N; j++ ) {
			if( frequentDrop[ j ] ) {
				distanceFrequentNodes[ i ][ j ] = distanceFromSourceL[ j ]; 
			} 
		}
		//cout<<distanceFrequentNodes[i].size()<<endl;
	}

	cout<<"SIZE OF FREQUENT NODES"<<distanceFrequentNodes.size()<<endl;
}

void queryInit(char *txtName) {
	degName = (char*)malloc(1+strlen(txtName) + 50);
	sprintf(degName, "%s.deg", txtName);
	
	inName = (char*)malloc(1+strlen(txtName) + 50);
	sprintf(inName, "%s.labelin", txtName);
	outName = (char*)malloc(1+strlen(txtName) + 50);
	sprintf(outName, "%s.labelout", txtName);

	timer tm;
	loadIndex();
	printf("load time %lf (ms)\n", tm.getTime()*1000); fflush(stdout);
}

/* returns the destination list given a list of passengers on the cab and a new passenger (source, destination) */
pair<vector< event >, double> midPointCalculation(int cabNode, vector< passenger > &passengerOnCab, passenger newPassenger) {
	int v = newPassenger.source;
	int w = newPassenger.destination;
	double cabv = query(cabNode, v); // should be 0 under cabNode == v
	vector< event > midPointList;
	double minDistanceTotal = MAX_DISTANCE;

	// Create eventIds
	vector< event > eventIds;
	for (int i=0; i<passengerOnCab.size(); i++) {
		event e(passengerOnCab[i].id, passengerOnCab[i].destination, true);
		eventIds.push_back(e);
	}
	event e(newPassenger.id, w, true);
	eventIds.push_back(e);

	// Try all permutation
	if (passengerOnCab.size() == 1) {
		int d0 = passengerOnCab[ 0 ].destination;
		double atd0 = passengerOnCab[ 0 ].actualTravelDist;
		double spd0 = passengerOnCab[ 0 ].shortestPathDist;

		double tvd0 = query(v, d0);	
		double td0w = query(d0, w);
		double tvw = query(v, w);
		double distanceTotal = 0;
		vector< event > candMidPointList;

		if ( ( atd0 + cabv + tvd0 <= alpha * spd0 ) && ( tvd0 + td0w <= alpha * tvw) ) {
			candMidPointList.clear();
			if (cabNode != v){  
				event e(newPassenger.id, v, false);
				candMidPointList.push_back( e ); 
			} 
			candMidPointList.push_back(eventIds[0]); candMidPointList.push_back(eventIds[1]);
			distanceTotal = cabv + tvd0 + td0w;
			if (distanceTotal < minDistanceTotal) {
				midPointList = candMidPointList;
				minDistanceTotal = distanceTotal;
			}
			
		}

		double twd0 = query(w, d0);
		if( atd0 + cabv + tvw + twd0 <= alpha * spd0) {
			candMidPointList.clear();
			if (cabNode != v) {
				event e(newPassenger.id, v, false);
				candMidPointList.push_back( e ); 
			}
			candMidPointList.push_back(eventIds[1]); candMidPointList.push_back(eventIds[0]);
			distanceTotal = cabv + tvw + twd0;
			if (distanceTotal < minDistanceTotal) {
				midPointList = candMidPointList;
				minDistanceTotal = distanceTotal;
			}
		}	

	}
	else if (passengerOnCab.size() == 2) {
		int s0 = passengerOnCab[ 0 ].source;
		int s1 = passengerOnCab[ 1 ].source;
		int d0 = passengerOnCab[ 0 ].destination;
		int d1 = passengerOnCab[ 1 ].destination;
	
		double ts0d0 = passengerOnCab[ 0 ].shortestPathDist;
		double ts1d1 = passengerOnCab[ 1 ].shortestPathDist;
		double atd0 = passengerOnCab[ 0 ].actualTravelDist;
		double atd1 = passengerOnCab[ 1 ].actualTravelDist;
		double tvd0 = query(v, d0);
		double tvd1 = query(v, d1);
		double td0d1 = query(d0, d1);
		double td1d0 = query(d1, d0);
		int secondPassengerFound = 0;

		double tvw = query(v, w);
		double td0w = query(d0, w);
		double twd0 = query(w, d0);
		double td1w = query(d1, w);
		double twd1 = query(w, d1);

		double distanceTotal = 0;
		vector< event > candMidPointList;

		/* 	1) V -> D First
			2) V -> D2 First  
			3) V -> W First */

		if( ( atd0 + cabv + tvd0 <= alpha*ts0d0) && ( atd1 + cabv + tvd0 + td0d1 <= alpha*ts1d1) && ( tvd0 + td0d1 + td1w <= alpha*tvw )  ) {
			candMidPointList.clear();
			if (cabNode != v) {
				event e(newPassenger.id, v, false);
				candMidPointList.push_back( e ); 
			}
			candMidPointList.push_back(eventIds[0]); candMidPointList.push_back(eventIds[1]); candMidPointList.push_back(eventIds[2]);
			distanceTotal = cabv + tvd0 + td0d1 + td1w;
			if (distanceTotal < minDistanceTotal) {
				midPointList = candMidPointList;
				minDistanceTotal = distanceTotal;
			}
		}
		
		if( ( atd0 + cabv + tvd0 <= alpha*ts0d0) && (tvd0 + td0w <= alpha*tvw) &&   (atd1 + cabv + tvd0 + td0w + twd1 <= alpha*ts1d1 ) ) {
			candMidPointList.clear();
			if (cabNode != v) {
				event e(newPassenger.id, v, false);
				candMidPointList.push_back( e ); 
			}
			candMidPointList.push_back(eventIds[0]); candMidPointList.push_back(eventIds[2]); candMidPointList.push_back(eventIds[1]);
			distanceTotal = cabv + tvd0 + td0w + twd1;
			if (distanceTotal < minDistanceTotal) {
				midPointList = candMidPointList;
				minDistanceTotal = distanceTotal;
			}
		}
		
		if ( (atd1 + cabv + tvd1 <= alpha*ts1d1) && (atd0 + cabv + tvd1 + td1d0 <= alpha*ts0d0) && ( tvd1 + td1d0 + td0w <= alpha*tvw ) ) {
			candMidPointList.clear();
			if (cabNode != v) {
				event e(newPassenger.id, v, false);
				candMidPointList.push_back( e ); 
			}
			candMidPointList.push_back(eventIds[1]); candMidPointList.push_back(eventIds[0]); candMidPointList.push_back(eventIds[2]);
			distanceTotal = cabv + tvd1 + td1d0 + td0w;
			if (distanceTotal < minDistanceTotal) {
				midPointList = candMidPointList;
				minDistanceTotal = distanceTotal;
			}
		}
		
		if ( (atd1 + cabv + tvd1 <= alpha*ts1d1) && (tvd1 + td1w <= alpha*tvw) &&  ( atd0 + cabv + tvd1 + td1w +twd0  <= alpha*ts0d0 ) ) {
			candMidPointList.clear();
			if (cabNode != v) {
				event e(newPassenger.id, v, false);
				candMidPointList.push_back( e ); 
			}
			candMidPointList.push_back(eventIds[1]); candMidPointList.push_back(eventIds[2]); candMidPointList.push_back(eventIds[0]);
			distanceTotal = cabv + tvd1 + td1w + twd0;
			if (distanceTotal < minDistanceTotal) {
				midPointList = candMidPointList;
				minDistanceTotal = distanceTotal;
			}
		}
		
		if( ( atd0 + cabv + tvw + twd0 <= alpha*ts0d0 ) && ( atd1 + cabv + tvw + twd0 + td0d1 <= alpha*ts1d1 ) ) {
			candMidPointList.clear();
			if (cabNode != v) {
				event e(newPassenger.id, v, false);
				candMidPointList.push_back( e ); 
			}
			candMidPointList.push_back(eventIds[2]); candMidPointList.push_back(eventIds[0]); candMidPointList.push_back(eventIds[1]);
			distanceTotal = cabv + tvw + twd0 + td0d1;
			if (distanceTotal < minDistanceTotal) {
				midPointList = candMidPointList;
				minDistanceTotal = distanceTotal;
			}
		}
		
		if( ( atd1 + cabv + tvw + twd1 <= alpha*ts1d1) && (atd0 + cabv + tvw + twd1 +td1d0 ) <= alpha*ts0d0 ) {
			candMidPointList.clear();
			if (cabNode != v) {
				event e(newPassenger.id, v, false);
				candMidPointList.push_back( e ); 
			}
			candMidPointList.push_back(eventIds[2]); candMidPointList.push_back(eventIds[1]); candMidPointList.push_back(eventIds[0]);
			distanceTotal = cabv + tvw + twd1 + td1d0;
			if (distanceTotal < minDistanceTotal) {
				midPointList = candMidPointList;
				minDistanceTotal = distanceTotal;
			}
		}
	}

	return make_pair(midPointList, minDistanceTotal);	
}

/* returns assignment of cabs to locations using Hungarian algorithm */
vector< int > assignFreeCabs( vector< int > &cabsLocation, vector< vector< passenger > > &passengerOnBoard,vector< passenger > &pickups, 
	int curTimeSlot ) {
	vector< vector<double> > costMatrix( cabsLocation.size() );
	
	cout<<" WE HAVE STARTED ASSIGNING"<<endl;
	cout<< costMatrix.size()<<" x "<<pickups.size()<<endl;
	for (int x = 0; x < costMatrix.size(); x++) {
		costMatrix[x].resize(pickups.size());
		for( int i = 0; i < pickups.size(); i++) {
			costMatrix[ x ][ i ] = MAX_DISTANCE;

			int v = pickups[ i ].source;
			int w = pickups[ i ].destination;
			
			if( !passengerOnBoard[ x ].size() ) {
				double dist = query(cabsLocation[ x ], v );
				if ( (curTimeSlot - pickups[i].startTimeSlot)*DELTA_TIME + dist / SPEED <= MINUTE_BEFORE_REJECTION ) {
					costMatrix[ x ][ i ] = dist;
				}
			}
			else {
				vector< event > midPointList = midPointCalculation(cabsLocation[ x ], passengerOnBoard[ x ], pickups[ i ]).first;

				if (midPointList.size() > 0) {
					double dist = query(cabsLocation[ x ], v );
					if ( (curTimeSlot - pickups[i].startTimeSlot)*DELTA_TIME + dist / SPEED <= MINUTE_BEFORE_REJECTION ) {
						costMatrix[ x ][ i ] = dist;
					}
				}
			}
		}
	}
	
	HungarianAlgorithm HungAlgo;
	vector<int> assignment;
	cout<<"DEBUG STATEMENT FOR HUNGARIAN ALGORITHM"<<endl;
	double cost = HungAlgo.Solve(costMatrix, assignment);
	cout<<"HUNGARIAN ALGORITHM RAN SUCCESSFULLY"<<endl;
	for( int i = 0; i < assignment.size(); i++) {
		if (assignment[ i ] != -1 ) {
			printf(" %d -> %d (%.4f)\n", i, assignment[i],costMatrix[ i ][ assignment[ i ] ]);
			if( costMatrix[ i ][ assignment[ i ] ] == MAX_DISTANCE ) {
				assignment[ i ] = -1;
			}
		}
	}
	return assignment;
}  

// Subset functions
// reduces the edges vector into edges only catering to nodes whose Lattitude < 40.6 (the marked ones)
// Also now the edges contain the new indices, which have nothing to do with original IdtoNode -> we need 'marked' for that conversion

void subset_edge(map<int, int> &marked, vector< vector<long long int> > &edges, vector< vector<double> > &edgeWeight) 
{
	vector< vector<long long int> > edgesTrunc;
	vector< vector<double> > edgeWeightTrunc;
	edgesTrunc.resize(marked.size());
	edgeWeightTrunc.resize(marked.size());

	for (int i=0; i<edges.size(); i++) 
	{
		int u = i;
		if (marked.find(u) == marked.end()) continue;
		
		for (int j=0; j<edges[i].size(); j++) 
		{
			int v = edges[i][j];
			double w = edgeWeight[i][j];
			if (marked.find(v) == marked.end()) continue;

			int mu = marked[u];
			int mv = marked[v];
			edgesTrunc[mu].push_back(mv);
			edgeWeightTrunc[mu].push_back(w);			
		}

	}
	edges.clear();
	edgeWeight.clear();
	edges = edgesTrunc;
	edgeWeight = edgeWeightTrunc;
}

void subset_dataset(map<int, int> &marked, vector< unordered_map<int, vector<request> > > &dataset) {
	vector< unordered_map<int, vector<request> > > datasetTrunc;
	datasetTrunc.resize(marked.size());
	for (int i=0; i<dataset.size(); i++) {
		int u = i;
		if (marked.find(u) == marked.end()) continue;
		for (unordered_map<int, vector<request> >::iterator iter=dataset[i].begin(); iter!=dataset[i].end(); iter++) {
			int timeSlot = iter->first;
			vector<request> requestList = iter->second;
			for (int j=0; j<requestList.size(); j++) {
				request req = requestList[j];
				int v = req.destination;
				if (marked.find(v) == marked.end()) continue;
				int mu = marked[u];
				int mv = marked[v];
				datasetTrunc[mu][timeSlot].push_back(req);
			}
			
		}
	}
	dataset.clear();
	dataset = datasetTrunc;
}


// Score calculation functions
void get_all_trips(int source,int start_time,vector<int > &tripList,
	vector< unordered_map<int, vector<request> > > &dataset ) {

	if( source < dataset.size() ) {
		if( dataset[source].find(start_time) != dataset[source].end() ) {

			for(int i = 0; i < dataset[ source][start_time].size(); i++) {
				tripList.push_back( dataset[source][ start_time ][ i ].destination );
			}

		}
	}
}

void get_all_trips_revenue(int source,int start_time,vector<pair<int,double> > &tripDestRev,
	vector< unordered_map<int, vector<request> > > &dataset ) {

	if( source < dataset.size() ) {
		if( dataset[source].find(start_time) != dataset[source].end() ) {

			for(int i = 0; i < dataset[ source][start_time].size(); i++) {
				tripDestRev.push_back( make_pair(dataset[source][ start_time ][ i ].destination, dataset[source][ start_time ][ i ].revenue) );
			}

		}
	}
}

void get_all_trips_fast(int source,int start_time,vector<request> &trips,
	vector< unordered_map<int, vector<request> > > &dataset ) {

	if( source < dataset.size() ) {
		if( dataset[source].find(start_time) != dataset[source].end() ) {
			trips = dataset[source][ start_time ];
		}
	}
}

double get_expected_trips(int v, vector<int> &d, int start_time, double alpha, vector<double> &distancesTravelled, double newDist, vector< double > &shortestPathDistances,
	vector< unordered_map<int, vector<request> > > &dataset, const vector<vector< double > > &distanceToDestination, const vector<vector< double > > &distanceFromDestination) {  	

	double expConst = 1.05;

	vector<request> tripList;
	get_all_trips_fast(v, start_time, tripList, dataset);

	if (tripList.size() == 0) {
		return 0;
	}

	// triplist has all the destinations
	int counter = 0;

	// calculating compatibility distance for the nodes in the dataset; counter returns it
	if (d.size() == 1) {
		double tvd, tsv;
		tsv = distancesTravelled[ 0 ] + newDist;
		tvd = distanceToDestination[ 0 ][ v ];

		for(int i = 0; i < tripList.size(); i++) {

			int w = tripList[ i ].destination;
			double tvw = query(v, w);
			double tdw = distanceFromDestination[ 0 ][ w ];

			if( tvd + tdw <= alpha*tvw ) {
				counter += 1; continue;
			}
			
			double twd = distanceToDestination[ 0 ][ w ];
			double distanceA = ( tsv + tvw + twd );
			double distanceB = alpha*shortestPathDistances[ 0 ];
			if( distanceA <= distanceB) {
				counter += 1; continue;
			}
		}

	}
	
	if (d.size() == 2) {

		double ts0v = distancesTravelled[ 0 ] + newDist;
		double ts1v = distancesTravelled[ 1 ] + newDist;

		double tvd0 = distanceToDestination[ 0 ][ v ];
		double tvd1 = distanceToDestination[ 1 ][ v ];

		double ts0d0 = shortestPathDistances[ 0 ];
		double ts1d1 = shortestPathDistances[ 1 ];

		double td0d1 = distanceFromDestination[ 0 ][ d[1] ];
		double td1d0 = distanceFromDestination[ 1 ][ d[0] ];

		for(int i = 0; i < tripList.size(); i++) {
			int w = tripList[ i ].destination;		

			double tvw = query(v, w);
			double td0w = distanceFromDestination[ 0 ][ w ];
			double twd0 = distanceToDestination[ 0 ][ w ];
			double td1w = distanceFromDestination[ 1 ][ w ];
			double twd1 = distanceToDestination[ 1 ][ w ];

			if( ( ts0v + tvd0 <= alpha*ts0d0) && ( ts1v + tvd0 + td0d1 <= alpha*ts1d1) && ( tvd0 + td0d1 + td1w <= alpha*tvw )  ) {
				counter += 1; continue;
			}
			if( ( ts0v + tvd0 <= alpha*ts0d0) && (tvd0 + td0w <= alpha*tvw) &&   (ts1v + tvd0 + td0w + twd1 <= alpha*ts1d1 ) ) {
				counter += 1; continue;
			}
			if ( (ts1v + tvd1 <= alpha*ts1d1) && (ts0v + tvd1 + td1d0 <= alpha*ts0d0) && ( tvd1 + td1d0 + td0w <= alpha*tvw ) ) {
				counter += 1; continue;
			}
			if ( (ts1v + tvd1 <= alpha*ts1d1) && (tvd1 + td1w <= alpha*tvw) &&  ( ts0v + tvd1 + td1w +twd0  <= alpha*ts0d0 ) ) {
				counter += 1; continue;
			}
			if( ( ts0v + tvw + twd0 <= alpha*ts0d0 ) && ( ts1v + tvw + twd0 + td0d1 <= alpha*ts1d1 ) ) {
				counter += 1; continue;
			}
			if( ( ts1v + tvw + twd1 <= alpha*ts1d1) && ( ts0v + tvw + twd1 +td1d0 ) <= alpha*ts0d0 ) {
				counter += 1; continue;
			}
		}
	}

	// printf("  %.4f - %.4f * %d = %.4f\n",newDist, pow(expConst, -newDist), counter, pow(expConst, -newDist) * counter);
	if(USE_DECAY)
		return pow(expConst, -newDist) * counter;
	else
		return counter;
	// return counter;
}

double get_scores_revenue(int v, int d, double distanceTravelled,
	int start_time, double alpha, vector< unordered_map<int, vector<request> > > &dataset,
	const vector< double > &distanceFromSource, const vector< double > &distanceToDestination, const vector< double > &distanceFromDestination) {  	

	vector<pair<int, double> > tripDestRev;
	get_all_trips_revenue(v, start_time, tripDestRev, dataset);

	double tsd, tvd, tsv;
	// no need ot calculate them for every iteration 
	if( tripDestRev.size() ) {
		tsd = distanceFromSource[ d ];
		tsv = distanceFromSource[ v ];
		tvd = distanceToDestination[ v ];
	}

	// calculating compatibility distance for the nodes in the dataset; counter returns it
	double sigmScore = 0;
	double sigmoidConst = 0.1;

	for(int i = 0; i < tripDestRev.size(); i++) 
	{
		int w = tripDestRev[ i ].first;
		double revenue = tripDestRev[ i ].second;

		// w = destination of the cab at midpoint  
		double tvw = query(v, w);
		double tdw = distanceFromDestination[ w ];
		double deltaDist = 1e10;

		// Original route: distanceTravelled + tvd
		// Route 1: distanceTravelled + tvd + tdw
		// Route 2: distanceTravelled + tvw + twd
		if( tvd + tdw <= alpha*tvw ) {
			deltaDist = min(deltaDist, tdw);
		}
		// cout<<"dist trav= "<<distanceTravelled<<endl;

		double twd = distanceToDestination[ w ];
		double distanceA = ( distanceTravelled + tvw + twd);
		double distanceB = alpha*tsd;
		if( distanceA <= distanceB ) {
			deltaDist = min(deltaDist, tvw + twd - tvd);
		}

		// This trip is not feasible
		if (deltaDist > 1e9) continue;

		// double profit = ( tsd + tvw ) - FUEL * ( distanceTravelled + tvd + tdw );
		double profit = revenue - FUEL * deltaDist;
		// cout<<"profit= "<<profit<<endl;
		// use direct e^x instead of that


		sigmScore += exp(sigmoidConst * profit);

		// cout<<"sigmoidConst= "<<sigmoidConst<<endl;
		// cout<<"sigmScore= "<< exp(sigmoidConst * profit) <<endl;
		// cout<<endl;
	}
	
	if(sigmScore==0)
		return 0;

	sigmScore= log(1+ sigmScore);
	// cout<<"final score= "<<sigmScore<<endl;

	// cout<<"midstop= "<<midStop<<endl;
	return sigmScore;
}

// ----> ith node to all nearby nodes and return their original node number along with the corresponding distance
vector<vector<pair<int, double> > > nearbyNodes;

void generateNearbyNodes(int N, const vector< vector<long long int> > &edges, const vector< vector<double> > &edgeWeight) { 

	vector<map<int, double> > checker(N);
	nearbyNodes.resize(N);

	int cntNearby = 0;

	// int visited[ dag.size() ];
	double distancesFromNode[ N ];
	for(int i=0;i<N;i++)
	{
		int node = i;
		for(int j=0;j<N;j++)
		{
			distancesFromNode[j]=MAX_DISTANCE;
		}
		distancesFromNode[ node ]=0;
		priority_queue<pair<double, int> > pq;
		pq.push(make_pair(0, node));
		// visited[ nodeToDagIndex[ node ] ] =1;

		while( !pq.empty() ) 
		{ 
			pair<double, int> x = pq.top(); 
			pq.pop();
			int u = x.second;

			if( distancesFromNode[ u ] > DELTA_DETOUR )
				break;

			if(checker[ node ].find(u) == checker[i].end() )
			{
				checker[ node ][ u ] = distancesFromNode[ u ];
				nearbyNodes[ node ].push_back(make_pair(u , distancesFromNode[ u ]));
				cntNearby++;
			}

			for(int j = 0; j < edges[ u ].size(); j++) { 	
				int v = edges[ u ][ j ];
				double alt = distancesFromNode[ u ] + edgeWeight[ u ][ j ];
				if( alt < distancesFromNode[ v ] )
				{ 	
					distancesFromNode[ v ] = alt;
					pq.push( make_pair( -alt, v) );
				}
			}
		}

	}

	printf("Finish generating nearby nodes (%d)\n",cntNearby); fflush(stdout);
}

double getClusterScore(int v, vector<int> &destinations, int start_time, double defaultAlpha, vector<double> &distancesTravelled, double newDist, vector< double > &shortestPathDistances,
	vector< unordered_map<int, vector<request> > > &dataset, const vector<vector< double > > &distanceToDestination, const vector<vector< double> > &distanceFromDestination) {	

	double scoreReturn=0;

	for(int i=0;i<=nearbyNodes[ v ].size();i++)
	{
		
		if(destinations.size()==1)
		{
			int u;
			double uDist;
			vector<request> trips;

			if(i == nearbyNodes[v].size())
			{
				u= v;
				uDist=0;
				get_all_trips_fast(u, start_time, trips, dataset);
			}
			else
			{
				pair<int, double> p = nearbyNodes[ v ][ i ];
				u = p.first;
				// u is the middle point
				double uDist = p.second;
				double timeOffset=0;

				if(TIME_VARYING)
					timeOffset = (uDist/SPEED)/MINUTE_PER_HISTORY_SLOT;	

				get_all_trips_fast(u, start_time + timeOffset, trips, dataset);

			}

			double shortestPath = shortestPathDistances[0];
			double allowed= defaultAlpha * shortestPath;

			double distTillU = distancesTravelled[ 0 ] + newDist + uDist;
			double sigmScore=0;
			double sigmoidConst = 0.1;
			double origDist=0;
			origDist = distanceToDestination[ 0 ][ v ];
			for(int j = 0; j < trips.size(); j++) 
			{
				int dest2 = trips[ j ].destination;
				double revenue = trips[ j ].revenue;

				// dest2 = destination of the cab at midpoint  
				double uToDest2 = query(u, dest2);
				double dest2ToDest1 = distanceToDestination[ 0 ][ dest2 ];
				double dest1ToDest2 = distanceFromDestination[ 0 ][ dest2 ];
				double deltaDist = 1e10;


				if( uToDest2 + dest2ToDest1 + distTillU <= allowed  ) {
					deltaDist = min(deltaDist, uDist+ uToDest2 + dest2ToDest1);

				}

				double uToDest1 = distanceToDestination[ 0 ][ u ];
				double distanceA = ( distTillU + uToDest1 + dest1ToDest2);
				double distanceB = defaultAlpha * uToDest2;
				if( distanceA <= distanceB && distanceA <= allowed) {
					deltaDist = min(deltaDist, uDist + uToDest1 + dest1ToDest2);
				}

				// This trip is not feasible
				if (deltaDist > 1e9) continue;

				scoreReturn += 1;
				// double profit = ( tsd + uToDest1 ) - FUEL * ( distanceTravelled + uToDest2 + dest1ToDest2 );
				//double profit = revenue - FUEL * (deltaDist-origDist);
				// use direct e^x instead of that
				//sigmScore += exp(sigmoidConst * profit);
			}
			// cout<<"sigmScore= "<<sigmScore<<endl;
			//scoreReturn+= exp(-uDist) * sigmScore;

		}

		if(destinations.size() == 2)
		{
			
			int u;
			double uDist;
			vector<request> trips;
			// cout<<"arre"<<endl;

			if(i == nearbyNodes[v].size())
			{
				u= v;
				uDist=0;  
				get_all_trips_fast(u, start_time, trips, dataset);
			}
			else
			{
				pair<int, double> p = nearbyNodes[ v ][ i ];
				u = p.first;
				// u is the middle point
				uDist = p.second;	
				double timeOffset=0;

				if(TIME_VARYING)
					timeOffset = (uDist/SPEED)/MINUTE_PER_HISTORY_SLOT;	

				get_all_trips_fast(u, start_time + timeOffset, trips, dataset);

			}

			int d1 = destinations[0];
			int d2 = destinations[1];
			double shortestPath1 = shortestPathDistances[0];
			double shortestPath2 = shortestPathDistances[1];

			// double allowed= defaultAlpha * shortestPathDistances;

			double distTillU1 = distancesTravelled[ 0 ] + newDist + uDist;
			double distTillU2 = distancesTravelled[ 1 ] + newDist + uDist;

			double sigmScore=0;
			double sigmoidConst = 0.1;			

			for(int j = 0; j < trips.size(); j++) 
			{
				int d3 = trips[ j ].destination;
				double revenue = trips[ j ].revenue;

				// dest2 = destination of the cab at midpoint  
				double shortestPath3 = query(u, d3);

				double dest3ToDest1 = distanceToDestination[ 0 ][ d3 ];
				double dest1ToDest3 = distanceFromDestination[ 0 ][ d3 ];
				double dest3ToDest2 = distanceToDestination[ 1 ][ d3 ];
				double dest2ToDest3 = distanceFromDestination[ 1 ][ d3 ];
				double dest2ToDest1 = distanceToDestination[ 0 ][ d2 ];
				double dest1ToDest2 = distanceFromDestination[ 0 ][ d2 ];

				double uToDest1 = distanceToDestination[ 0 ][ u ];
				double uToDest2 = distanceFromDestination[ 1 ][ u ];
				double uToDest3 = shortestPath3;
				
				double origDist=0;
				origDist = distanceToDestination[ 0 ][ v ] + dest1ToDest2;

				double deltaDist = 1e10;

				double allowed1 = defaultAlpha * shortestPath1;
				double allowed2 = defaultAlpha * shortestPath2;
				double allowed3 = defaultAlpha * shortestPath3;
				// see for permutations


				// 123
				if((distTillU1 + uToDest1 <=allowed1) && (distTillU2 + uToDest1 + dest1ToDest2 <=allowed2) && (uToDest1 + dest1ToDest2 + dest2ToDest3 <= allowed3))
				{
					deltaDist = min(deltaDist, uDist + uToDest1 + dest1ToDest2 + dest2ToDest3);
				}
				// 132
				if((distTillU1 + uToDest1 <=allowed1) && (distTillU2 + uToDest1 + dest1ToDest3 + dest3ToDest2 <=allowed2) && (uToDest1 + dest1ToDest3 <= allowed3))
				{
					deltaDist = min(deltaDist, uDist + uToDest1 + dest1ToDest3 + dest3ToDest2);
				}
				// 213
				if((distTillU1 + uToDest2 + dest2ToDest1 <=allowed1) && (distTillU2 + uToDest2 <=allowed2) && (uToDest2 + dest2ToDest1 + dest1ToDest3 <= allowed3))
				{
					deltaDist = min(deltaDist, uDist + uToDest2 + dest2ToDest1 + dest1ToDest3);
				}
				// 231
				if((distTillU1 + uToDest2 + dest2ToDest3 + dest3ToDest1 <=allowed1) && (distTillU2 + uToDest2 <=allowed2) && (uToDest2 + dest2ToDest3 <= allowed3))
				{
					deltaDist = min(deltaDist, uDist + uToDest2 + dest2ToDest3 + dest3ToDest1);
				}
				// 312
				if((distTillU1 + uToDest3 + dest3ToDest1 <=allowed1) && (distTillU2 + uToDest3 + dest3ToDest1 + dest1ToDest2 <=allowed2) && (uToDest3 <= allowed3))
				{
					deltaDist = min(deltaDist, uDist + uToDest3 + dest3ToDest1 + dest1ToDest2);
				}
				// 321
				if((distTillU1 + uToDest3 + dest3ToDest2 + dest2ToDest1 <=allowed1) && (distTillU2 + uToDest3 + dest3ToDest2 <=allowed2) && (uToDest3 <= allowed3))
				{
					deltaDist = min(deltaDist, uDist + uToDest3 + dest3ToDest2 + dest2ToDest1);
				}



				// This trip is not feasible
				if (deltaDist > 1e9) continue;

				scoreReturn += 1;
				// double profit = ( tsd + uToDest1 ) - FUEL * ( distanceTravelled + uToD3 + dest1ToDest2 );
				//double profit = revenue - FUEL * (deltaDist-origDist);
				// use direct e^x instead of that
				//sigmScore += exp(sigmoidConst * profit);
			}
			//scoreReturn+= exp(-uDist) * sigmScore;


		}

	}	

	if(scoreReturn==0) {
		// printf("Get Cluster Score at %d = %.2f\n",v, scoreReturn); fflush(stdout);
		return 0;
	}

	// printf("Before log %d = %.2f... ",v, scoreReturn); fflush(stdout);
	//scoreReturn= log(1+ scoreReturn);

	// printf("Get Cluster Score at %d = %.2f\n",v, scoreReturn); fflush(stdout);
	return scoreReturn;
	
}


// actualTravelDist = actual from source to midStop
int get_expected_trips2(int source, pair< int, int> destination, int midStop,
	int start_time, double alpha, vector<double> varAlphas, vector< unordered_map<int, vector<request> > > &dataset,
	const vector< double > &distanceFromSource, const vector< double > &distanceToDestination, const vector< double > &distanceFromDestination, double actualTravelDist = 0.0) {  	

	vector<int> tripList;
	
	get_all_trips(midStop, start_time, tripList, dataset);

	int counter = 0;
	double alpha1 = varAlphas[1]; double alpha2 = varAlphas[2];
	int s = source;
	int d = destination.first;
	int d2 = destination.second;
	int v = midStop;
	double tvd, tsv, ts_d2, td_d2, td2_d, tsd, tv_d2;
	// no need to calculate them for every iteration
	if( tripList.size() ) {
		tsd = distanceFromSource[ d ];
		tsv = distanceFromSource[ v ];
		ts_d2 = distanceFromSource[ d2 ];
		tvd = distanceToDestination[ v ];
		td_d2 = distanceFromDestination[ d2 ];
		td2_d = distanceToDestination[ d2 ];
		tv_d2 = query(v, d2);
	}

	for(int i = 0; i < tripList.size(); i++) {

		int w = tripList[ i ];
		double tvw = query(v, w);
		double td2_w = query(d2, w);
		double tw_d2 = query(w, d2);
		double tdw = distanceFromDestination[ w ];
		double twd = distanceToDestination[ w ];

		/* 	1) V -> D First
			2) V -> D2 First  
			3) V -> W First */
		if(  ( actualTravelDist + tvd + td_d2 <= alpha2*ts_d2 && ( tvd + td_d2 + td2_w <= alpha*tvw ) ) 
					|| (tvd + tdw <= alpha*tvw) &&  (actualTravelDist + tvd + tdw +tw_d2 ) <= alpha2*ts_d2  ) {
			counter += 1;
		}
		else if( ( actualTravelDist + tv_d2 + td2_d <= alpha1*tsd && ( tv_d2 + td2_d + tdw <= alpha*tvw ) ) 
					|| ( (tv_d2 + td2_w <= alpha*tvw) &&  ( actualTravelDist + tv_d2 + td2_w +twd  <= alpha1*tsd ) ) ) {
			
			counter += 1;
		}
		else if( ( actualTravelDist + tvw + twd <= alpha1*tsd && ( actualTravelDist + tvw + twd + td_d2 <= alpha2*ts_d2 ) ) 
					|| ( actualTravelDist + tvw + tw_d2 <= alpha2*ts_d2) &&  ( actualTravelDist + tvw + tw_d2 +td2_d ) <= alpha1*tsd ) {
			
			counter += 1;
		}
	}

	return counter;
}

/*
Dijkstra algorithm
Given (S,D), return the shortest path
*/
vector<int> dijkstra_lengths(int N, int S, int D, vector< double > &distanceFromSource,
	const vector< vector<long long int> > &edges, const vector< vector<double> > &edgeWeight) { 

	vector<int> prevNode(N);
	for(int i = 0; i < N; i++)
	{ 	distanceFromSource[ i ] = MAX_DISTANCE;
		prevNode[ i ] = -1;
	}

	distanceFromSource[ S ] = 0;	
	priority_queue< pair<float, int> > dj;
	dj.push( make_pair(0, S) );
	
	pair<float, int> x;
	int u, v;
	float alt;

	while( dj.size() ) 
	{ 
		x = dj.top(); 
		dj.pop();
		u = x.second;

		if( distanceFromSource[ u ] >= MAX_DISTANCE )
			break;

		for(int i = 0; i < edges[ u ].size(); i++) { 	
			v = edges[ u ][ i ];
			alt = distanceFromSource[ u ] + edgeWeight[ u ][ i ];
			if( alt < distanceFromSource[ v ] )
			{ 	distanceFromSource[ v ] = alt;
				dj.push( make_pair( -alt, v) );
				prevNode[ v ] = u;
			}
		}
	}

	for(int i = 0; i < N; i++)
	{ 	if( distanceFromSource[ i ] >= MAX_DISTANCE ) {
			distanceFromSource[ i ] = MAX_DISTANCE;
			prevNode[ i ] = -1;
		}
	}

	vector<int> path;
	int node = D;
	while( true ) {
		path.push_back( node );
		if( ( node == S ) || ( prevNode[ node ] == -1) )
			break;
		node = prevNode[ node ];
	}
	
	reverse(path.begin(), path.end());

	return path;
}

/*
Greedy algorithm:
Given (S,D), take a step where node score / edge weight is maximized. Return a path
*/
vector<int> greedyScoreDivDist(int N, int S, int D, double alpha, vector< int > &expectedTrips,
	const vector< double > &distanceFromSource, vector< vector<long long int> > &edges, vector< vector<double> > &edgeWeight) { 

	vector<float> dist(N);
 
	for(int i=0; i< N ; i++) {
		dist[ i ] = 100000;
	}

	dist[S] = 0;

	// ( (score, node), (distance, path) )
	priority_queue< pair< pair<float, int>, pair< float, vector<int> > > > dj;
	
	vector<int> path;
	path.push_back(S);
	dj.push( make_pair( make_pair(0, S), make_pair(0, path) ) );

	while (dj.size()) { 
		pair< pair<float, int>, pair< float, vector<int> >  > x = dj.top();
		dj.pop();
		
		int u = x.first.second;
		//printf(" At %lld: %.4f [%.4f]\n",u,-x.first.first, x.second.first);
		
		if( x.second.first > alpha * distanceFromSource[D]  )
			continue; 
		
		if( u == D ) {
			return x.second.second;
		}
		
		for(int i=0; i< edges[u].size(); i++) {
			int v = edges[u][i];

			if (find(x.second.second.begin(), x.second.second.end(), v) != x.second.second.end()) {
				continue;
			}

			float alt = -x.first.first + ( edgeWeight[u][i] / (expectedTrips[v]+1) );
			if( (x.second.first + edgeWeight[u][i]) < dist[v] )
			{ 	
				dist[v] = (x.second.first + edgeWeight[u][i]); 
				path = x.second.second;
				path.push_back(v);
				dj.push( make_pair( make_pair(-alt, v), make_pair(dist[v], path) ) );
			}
		}
	}

	// not reachable
	path.clear();
	return path; 

}

pair<double, vector<int> > dfsUtil(int N, int source, int v, vector<int> &destinations, double defaultAlpha, vector<int> &recStack, int timeSlot, double allowed, vector<double> &actualTravelDists, vector<double> &shortestPathDistances, const vector< vector<long long int> > &edges,
 		const vector< vector<double> > &edgeWeight,  const vector< double > &distanceFromSource, const vector<vector< double > > &distanceToDestination, const vector<vector< double > > &distanceFromDestination, vector< unordered_map<int, vector<request> > > &dataset)
{
	int destination = destinations[0];
	double travelled = actualTravelDists[0];
	vector<int> bestPath;

	if(travelled + distanceToDestination[0][v] > allowed)
	{
		if(source==v)
		{
			cout<<"  trouble"<<endl<<endl<<endl;
			cout<<"distance= "<<distanceToDestination[0][v]<<endl;
			cout<<"travelled= "<<travelled<<endl;
		} 

		return make_pair(-1, bestPath);
	}
	double timeOffset= 0;
	if(TIME_VARYING)
		timeOffset= ((int)(travelled/SPEED))/MINUTE_PER_HISTORY_SLOT;

	double node_score = get_expected_trips(v, destinations, timeSlot + timeOffset, defaultAlpha, actualTravelDists, 0 , shortestPathDistances, 
						dataset, distanceToDestination, distanceFromDestination);
	if (source == v) 
		node_score = 0;
	
	if(v == destination) {
		return make_pair(node_score, bestPath);
	}

	double maxScoreFromNode = 0;
	int maxNode = -1;
	

	recStack[v]=1;


	for(int i=0;i< edges[v].size();i++)
	{
		int neighbour=  (int) edges[v][i];
		double edgeDist = edgeWeight[v][i];
		if(travelled + edgeDist <= allowed && recStack[neighbour] == 0 )
		{
			actualTravelDists[0]+=edgeDist;
			pair<double, vector<int> > result = dfsUtil(N, source, neighbour, destinations, defaultAlpha, recStack, timeSlot, allowed, actualTravelDists, shortestPathDistances, edges, 
							edgeWeight, distanceFromSource, distanceToDestination, distanceFromDestination, dataset );
			if (result.first > maxScoreFromNode) {
				maxScoreFromNode = result.first;
				bestPath = result.second;
				maxNode = neighbour;
			}
			actualTravelDists[0]-=edgeDist;
		}

	}

		
	// cout<<"v= "<<v<<" score here= "<<scoreHere<<endl;
	
	recStack[v]=0;

	bestPath.push_back(maxNode);
	return make_pair(node_score + maxScoreFromNode, bestPath);

}

pair<double, vector<int> > optimalScore(int N, int source, vector<int> destinations, int timeSlot, double defaultAlpha, vector<double> &actualTravelDists, vector<double> &shortestPathDistances, vector< unordered_map<int, vector<request> > > &dataset,
	const vector< vector<long long int> > &edges, const vector< vector<double> > &edgeWeight, const vector< double > &distanceFromSource, const vector<vector< double > > &distanceToDestination, const vector<vector< double > > &distanceFromDestination)
{
	double shortestPathLen = shortestPathDistances[0];
	double allowed = defaultAlpha * shortestPathLen ;
	// cout<<shortestPathLen<<endl;
	// cout<<"allowed "<<allowed<<endl;
	vector<int> recStack(N);
	recStack[source]=1;
	cout<<"source= "<<source<<endl;
	// cout<<edges.size()<<endl;
	int destination= destinations[0];
	cout<<edges[source].size()<<endl;
	cout<<actualTravelDists[0];
	pair<double, vector<int> > result = dfsUtil(N, source, source, destinations, defaultAlpha, recStack, timeSlot, allowed, actualTravelDists, shortestPathDistances , edges, 
							edgeWeight, distanceFromSource, distanceToDestination, distanceFromDestination, dataset);

	cout << "Optimal score: " << result.first << endl;
	vector<int> path = result.second;
	path.push_back(source);
	reverse(path.begin(), path.end());

	return make_pair(result.first, path);
}


/*
DAG score algorithm
Given (S,D), construct a DAG based on distance to destination, then perform a dynamic programming to obtain max score path
*/
pair<double, vector<int> > findDAGPath(int N, int source, vector<int> destinations, int timeSlot, double defaultAlpha, vector<double> &actualTravelDists, vector<double> &shortestPathDistances, vector< unordered_map<int, vector<request> > > &dataset,
	const vector< double > &distanceFromSource, const vector<vector< double > > &distanceToDestination, const vector<vector< double > > &distanceFromDestination, const vector< vector<long long int> > &edges, const vector< vector<double> > &edgeWeight) 
{
	int firstDestination = destinations[0];

	vector< pair< double, int> > dag;

	double remainShortestDist = 0;
	bool alphaImpossible = false;
	for(int i=0; i < N ; i++) {
		remainShortestDist = 0;
		alphaImpossible = false;
		for (int j=0; j<destinations.size(); j++) {
			if (j) remainShortestDist += distanceToDestination[j][destinations[j-1]];
			// if (i == destinations[j] || i == source) printf("%d: %.2f + %.2f + %.2f + %.2f > %.2f * %.2f ? %d\n",i,actualTravelDists[j],distanceFromSource[i],distanceToDestination[0][i],remainShortestDist,defaultAlpha,shortestPathDistances[j],actualTravelDists[j] + distanceFromSource[i] + distanceToDestination[0][i] + remainShortestDist > defaultAlpha * shortestPathDistances[j]); fflush(stdout);
			if (actualTravelDists[j] + distanceFromSource[i] + distanceToDestination[0][i] + remainShortestDist > defaultAlpha * shortestPathDistances[j]) {
				alphaImpossible = true; break;
			}
		}
		
		if (!alphaImpossible) dag.push_back( make_pair(-distanceToDestination[0][ i ], i ) );
	}

	printf("Done! %d\n",(int)dag.size()); fflush(stdout);

	vector<int> path;
	if(dag.size() <= 1)
	{
		// This happens because query is not as accurate as Dijkstra
		vector <double> dummy(N);
		path= dijkstra_lengths(N, source, firstDestination, dummy, edges, edgeWeight);
		cout<<"shortestpaths used instead DAG" <<endl;
		return make_pair(-1, path);
	}


	sort( dag.begin(), dag.end() );

	//here initialise all the distances (nearby nodes)


	// (Distance, Last node) at dag[i] given a particular score
	// ----------------------------------------------------------------------
	// vector< map<int, pair<double, pair<long long int, long long int> > > > scores ( dag.size() );
	vector< vector < pair< double, double > > > scores(dag.size());
	vector< vector <pair< int, int> > > store(dag.size());
	// --------------------------------------------------------------------
	
	vector<int> nodeToDagIndex(N);
	vector<int> dagIndexToNode(N);
	// cout<<dag.size()<<endl;
	for (int i=0; i<N; i++) {
		nodeToDagIndex[i] = -1;
		dagIndexToNode[i] = -1;
	}

	// for nodes who are in the dag after filtering -> new index
	for(int i=0; i < dag.size() ; i++) {
		nodeToDagIndex[ dag[i].second ] = i;
		dagIndexToNode[ i ]= dag[i].second;
	}
		
	int startIndex = nodeToDagIndex[source];
	int endIndex = nodeToDagIndex[firstDestination];

	remainShortestDist = 0;
	double delta_dist = 1e10;
	for (int j=0; j<destinations.size(); j++) {
		if (j) remainShortestDist += distanceToDestination[j][destinations[j-1]];
		delta_dist = min(delta_dist, defaultAlpha * shortestPathDistances[j] - actualTravelDists[j] - remainShortestDist);
	}
	// double delta_dist= defaultAlpha * shortestPathDistances[0] - actualTravelDists[0];

	// cout<<"delta dist= "<<delta_dist<<endl;

	// scores array is the dp
	// scores[ startIndex ][ 0 ] = make_pair(0, make_pair(-1,-1));
	double neginf = -1e10;
	int negint = INT_MIN;

	vector<pair<double, double> > forscore;
	vector<pair<int, int> > forstore;
	for(int j=0; j<=DISCRETE_TIME; j++)
	{
		forscore.push_back( make_pair(neginf, neginf ) );
		forstore.push_back( make_pair(negint, negint ));
	}

	for(int i=0; i<dag.size(); i++)
	{
		scores[i] = forscore;
		store[i] = forstore;
	}

	scores[startIndex][ 0 ].first=0;
	scores[startIndex][ 0 ].second=0;
	store[startIndex][0].first =-1;
	store[startIndex][0].second =-1;
	
		
	for(int i = startIndex; i < dag.size(); i++) 
	{
			
		for(int k=0;k<DISCRETE_TIME;k++)
		{

			int u = dag[i].second;
			if(scores[i][k].first<0)
				continue;

			for(int j = 0; j < edges[u].size(); j++) 
			{
				int v = edges[u][j];
				int vIndex = nodeToDagIndex[ v ];
				if( nodeToDagIndex[ v ] == -1 )
					continue;
				//some nodes are adjacent to themselves in the graph
				if(u == v)
					continue;

				double edge_dist= edgeWeight[ u ][ j ];
				double prev_dist= scores[i][k].second;
				// needs to change -> index wise


				// assuming dist in km and time in min and speed is 60 kmph
				// double edge_time= edge_dist;
				int newTime= DISCRETE_TIME * (edge_dist + prev_dist)/delta_dist;
				// int offset= (DISCRETE_TIME*edge_dist/delta_dist);

				if(newTime <= DISCRETE_TIME && vIndex > i )
				{

					double newDist = edge_dist + prev_dist;

					map< int, int > prevSource;
					int timeOffset=0;

					if(TIME_VARYING)
						timeOffset = ((int)(newDist/SPEED))/MINUTE_PER_HISTORY_SLOT;

					double node_score = 0;

					if(USE_CLUSTER)
					{
						node_score = getClusterScore(v, destinations, timeSlot + timeOffset, defaultAlpha, actualTravelDists, newDist, shortestPathDistances, 
						dataset, distanceToDestination, distanceFromDestination);
					}
					else
					{
						node_score = get_expected_trips(v, destinations, timeSlot + timeOffset, defaultAlpha, actualTravelDists, newDist, shortestPathDistances, 
						dataset, distanceToDestination, distanceFromDestination);
					}

					
					//cout<<"dag score= "<<"v= "<<v<<" "<<node_score<<endl;
					// double node_score = getClusterScore(v, destinations, timeSlot + ((int)(newDist/SPEED))/MINUTE_PER_HISTORY_SLOT, defaultAlpha, actualTravelDists, newDist, shortestPathDistances, 
						// dataset, distanceToDestination, distanceFromDestination);
					// if (node_score_old > 0) printf(" Node[%d]: prev %lld now %.4f\n", v, node_score_old, node_score);

					if((scores[ nodeToDagIndex[ v ] ][ newTime ].first  <  scores[ i ][ k ].first +  node_score) || 
						((scores[ nodeToDagIndex[ v ] ][ newTime ].first  ==  scores[ i ][ k ].first +  node_score) && newDist < scores[ nodeToDagIndex[ v ] ][ newTime ].second))
					{
						scores[ nodeToDagIndex[ v ] ][ newTime ].first = scores[ i ][ k ].first + node_score;
						scores[ nodeToDagIndex[ v ] ][ newTime ].second = newDist;
						store[ vIndex ][ newTime ]= make_pair(i, k);
					}				

				}
		
			
			}
		}
		
	}

	// for(int i=0; i<dag.size(); i++)
 //    {
 //        // vector <double> fogRow;
 //        for(int j=0; j<DISCRETE_TIME; j++)
 //        {
 //             if(scores[i][j]>=0)
 //             {
 //             	cout<<"i= "<<i<<" j= "<<j<<" sc= ";
 //             	cout<<scores[i][j]<<" ";
 //             }
 //        }
 //        // scores.push_back(fogRow);
 //        cout<<endl;
 //    }


	long long int bestScore = 0;


	// backtracking -> needed to change just like updation
	double max_score=-1;
	// long long int maxindex=-1;
	int max_time=-1;
	for( int popp=0;popp<=DISCRETE_TIME;popp++) 
	{
		//printf(" -- Score = %d (%.4f < %.4f?)\n",it->first, it->second.first, distFromSourceToDestination*alpha);
		if( scores[endIndex][popp].first>max_score ) 
		{
			max_score= scores[endIndex][popp].first;
			// maxindex=it->second.second.first;
			max_time= popp;
		}
	}

	if(max_score==-1)
	{
		vector <double> dummy(N);
		path= dijkstra_lengths(N, source, firstDestination, dummy, edges, edgeWeight);
		cout<<"shortestpaths used instead DAG" <<endl;
		return make_pair(-1, path);
	}

	int traceLocation = endIndex;
	int tracetime = max_time;

	//check for max_score -1
	// in that case should we return the shortest path

	while (traceLocation != -1) 
	{
		//printf(" At %lld[%lld]: Score = %d (dist: %.4f)\n",traceLocation, nodeToDagIndex[traceLocation], traceScore, scores[ nodeToDagIndex[traceLocation] ][traceScore].first);
		path.push_back(dagIndexToNode[ traceLocation] );
		int prevtime = store[ traceLocation ][ tracetime ].second;
		traceLocation = store[ traceLocation ][tracetime].first;
		tracetime = prevtime;
	}
	reverse(path.begin(), path.end());

	cout<<"BEST SCORE:   "<<max_score<<endl;
	double pathlen=getPathDist(path, edges, edgeWeight);
	cout<<"Pathlen= "<<pathlen<<endl;
	cout<<"delta_dist= "<<delta_dist<<endl;
	// for(int i=0;i<path.size();i++ )
	// {
	// 	cout<<path[i]<<endl;
	// }
	double pathlen2=distanceFromSource[firstDestination];
	cout<<"disjkstra Pathlen= "<<pathlen2<<endl;

	if (DEBUG && dag.size() <= 40) {
		printf("s:%d e:%d\n",startIndex, endIndex);
		for(int i = 0; i < dag.size(); i++) {
			int u = dag[i].second;
			double node_score = get_expected_trips(u, destinations, timeSlot, defaultAlpha, actualTravelDists, 0, shortestPathDistances, 
						dataset, distanceToDestination, distanceFromDestination);
			printf("%d %.4f\n",i,node_score);

		}
		for (int i = 0; i<dag.size(); i++) {
			int u = dag[i].second;
			for(int j = 0; j < edges[u].size(); j++) {
				int v = edges[u][j];
				if( nodeToDagIndex[ v ] == -1 )
					continue;
				//some nodes are adjacent to themselves in the graph
				if(u == v)
					continue;
				printf("%d %d %.4f\n",nodeToDagIndex[u], nodeToDagIndex[v],edgeWeight[ u ][ j ]);
			}
		}
		printf("Path: ");
		for (int i=0; i<path.size(); i++) {
			printf(" %d",nodeToDagIndex[path[i]]);
		}
		printf("\n");
	}

	return make_pair(max_score,path);

}

/*
DAG extended score algorithm
Given (S,D), construct a DAG based on distance to destination, then perform a dynamic programming to obtain max score path. 
Back edges are allowed as long as the distance is less than the maxDepth parameter.
Cycle is not allowed.
For multi passenger cab, we always seek for planning the best route to the first destination
*/

vector< vector< pair<double, vector<int> > > > extendEdge;

void extendEdges(int source, int node, double maxDepth, vector<int> &path, vector< pair<double, vector<int> > > &paths,
	double pathLen, const vector< vector<long long int> > &edges, const vector< vector<double> > &edgeWeight) {
	
	for( int i = 0; i < path.size(); i++) {
		if( node == path[i] )
			return;
	}
	path.push_back( node );

	if( (pathLen >= maxDepth) && ( path.size() > 2 ) ) {
		path.pop_back();
		return ;
	}

	if( path.size() > 1 ) {
		paths.push_back( make_pair(pathLen, path) ) ;
	}
	
	for(int j=0; j < edges[ node ].size(); j++) {
		int newNode = edges[ node ][j];
		extendEdges( source, newNode, maxDepth, path, paths, pathLen + edgeWeight[ node ][ j ], edges, edgeWeight ) ;
	}
	//cout<<"node= "<<node<<endl;

	path.pop_back();
}

 
void assignExtendEdge(int N, double maxDepth, const vector< vector<long long int> > &edges, const vector< vector<double> > &edgeWeight) {
	extendEdge.resize(N);
	int cnt_extended = 0;
	for( int i = 0; i < N; i++) {
		vector<int> path;
		vector< pair< double, vector<int> > > paths;
		extendEdges(i, i, maxDepth, path, paths, 0, edges, edgeWeight);
		cnt_extended += paths.size();
		extendEdge[ i ] = paths;
	}
	printf("Extend edge assigned (Total = %d)\n",cnt_extended); fflush(stdout);
}

int short_outputs=0;

pair<double, vector<int> >findDAGExtendedPath(int N, int source, vector<int> destinations, int timeSlot, double defaultAlpha, vector<double> &actualTravelDists, vector<double> &shortestPathDistances, double maxDepth, vector< unordered_map<int, vector<request> > > &dataset,
	const vector< double > &distanceFromSource, const vector<vector< double > > &distanceToDestination, const vector<vector< double > > &distanceFromDestination, const vector< vector<long long int> > &edges, const vector< vector<double> > &edgeWeight) {


	int firstDestination = destinations[0];

	map<int, int> dummy;

	vector< vector< long long int> > extendEdgeWeights(N);
	for( int i = 0; i < N; i++) {
		extendEdgeWeights[i].resize(extendEdge[i].size());
	}

	vector< pair< double, int> > dag;

	double remainShortestDist = 0;
	bool alphaImpossible = false;
	for(int i=0; i < N ; i++) {
		remainShortestDist = 0;
		alphaImpossible = false;
		for (int j=0; j<destinations.size(); j++) {
			if (j) remainShortestDist += distanceToDestination[j][destinations[j-1]];
			if (actualTravelDists[j] + distanceFromSource[i] + distanceToDestination[0][i] + remainShortestDist > defaultAlpha * shortestPathDistances[j]) {
				alphaImpossible = true; break;
			}
		}
		
		if (!alphaImpossible) dag.push_back( make_pair(-distanceToDestination[0][ i ], i ) );
	}

	printf("Done! %d\n",(int)dag.size()); fflush(stdout);

	vector<int> path;
	if(dag.size() <= 1)
	{
		// This happens because query is not as accurate as Dijkstra
		vector <double> dummy(N);
		path= dijkstra_lengths(N, source, firstDestination, dummy, edges, edgeWeight);
		cout<<"shortestpaths used instead DAG" <<endl;
		return make_pair(-1,path);
	}

	sort( dag.begin(), dag.end() );

	// (Distance, (Last Node, Ex Edge used) ) at dag[i] given a particular score
	vector< vector <pair<double, double> > > scores( dag.size() );
	vector< vector< pair<vector<int> ,int > > > store( dag.size() );

	
	vector< int> nodeToDagIndex(N);
	vector< int> dagIndexToNode(N);
	for (int i=0; i<N; i++) {
		nodeToDagIndex[i] = -1;
		dagIndexToNode[i] = -1;
	}
	for(int i=0; i <  dag.size() ; i++) {
		nodeToDagIndex[ dag[i].second ] = i;
		dagIndexToNode[ i ] = dag[i].second;
	}
	// cout<<"start 4"<<endl;
	// fflush(stdout);
	
	int startIndex = nodeToDagIndex[ source ];
	int endIndex = nodeToDagIndex[ firstDestination ];
	
	double neginf = -1e10;
	int negint = INT_MIN;

	remainShortestDist = 0;
	double delta_dist = 1e10;
	for (int j=0; j<destinations.size(); j++) {
		if (j) remainShortestDist += distanceToDestination[j][destinations[j-1]];
		delta_dist = min(delta_dist, defaultAlpha * shortestPathDistances[j] - actualTravelDists[j] - remainShortestDist);
	}

	vector<int> tadd;
	tadd.push_back(negint);
	vector <pair<double, double> > forRow;
	vector <pair<vector<int>, int > > forstore;

	for(int j=0; j<=DISCRETE_TIME; j++)
	{
		 forRow.push_back( make_pair(neginf, neginf) );
		 forstore.push_back( make_pair(tadd, negint ));
	}
	for(int i=0; i<dag.size(); i++)
	{       
		scores[i] = forRow;
		store[i] = forstore;
	}

	scores[startIndex][ 0 ]=make_pair(0, 0);
	store[startIndex][0].first =tadd;
	store[startIndex][0].second =-1;

	for(int i = startIndex; i < dag.size(); i++) 
	{
		for(int k=0;k<=DISCRETE_TIME;k++)
		{
			int u = dag[i].second;
			if(scores[i][k].first < 0)
				continue;
			// cout << i << " " << k << endl;
			for(int j = 0; j < extendEdge[u].size(); j++) 
			{
				
				pair<double, vector<int> > v = extendEdge[u][j];
				bool continueLoop = false;
				bool invalidNode = false;
				// bool reachAbleError = false;
				int pathSize = v.second.size();
				vector<int> ex_edge = v.second;

				if( pathSize < 1)
					continue;

				int lastNode =  ex_edge[ pathSize - 1 ];
				if(nodeToDagIndex[ lastNode ]<= i)
					continue;

				if ( lastNode == u )
					continue;

				// All nodes in the extended path except the last node has to be before u in the DAG
				for(int piter = 1; piter < pathSize ; piter++ ) 
				{
					int pathNode = ex_edge[ piter]; 
					
					if( ( nodeToDagIndex[ pathNode ] >= i )  ^  ( piter == ( pathSize -1 ) ) ) {
						continueLoop = true;
						break;
					}

					if( nodeToDagIndex[ pathNode ] == -1 ) {
						invalidNode = true;
						break;
					}
				}

				if( continueLoop || invalidNode )
					continue;

				map< int, int> markNodes;
				for(int k1 = 1; k1 < pathSize ; k1++ ) {
					markNodes[ ex_edge[ k1 ] ] = 1;
				} 

				int vIndex = nodeToDagIndex[ lastNode ];

				bool reachAbleError = false;
				int traceLoc = nodeToDagIndex[ u ];
				int trace_time = k;

				// bscktracking checking
				while ( (traceLoc != -1 ) && ( !reachAbleError ) ) 
				{
					// TODO: pass by reference
					pair<vector<int>, int > prev = store[ traceLoc ][ trace_time ];
					vector<int> prev_loc = prev.first;
					int prev_time = prev.second;
					// cout<<"huh4"<<endl;
					for(int loc = prev_loc.size()-1 ; loc>=0 ; loc-- )
					{
						if( loc >= 0 && markNodes.find( prev_loc[ loc ]  ) != markNodes.end() ) 
						{
							reachAbleError = true;
							break;
						}
					}

					// cout<<"traceLoc= "<<traceLoc<<endl;
					if (traceLoc != -1) {
						trace_time = prev_time; 
					}
					// cout<<"tracetime= "<<trace_time<<endl;
					if(prev_loc[0] < 0)
						break;
					else
						traceLoc = nodeToDagIndex [ prev_loc[0] ];

					// cout<<"traceLoc2= "<<traceLoc<<endl;
					if(traceLoc != -1)
						if( distanceToDestination[0][ prev_loc[ 0 ] ] - distanceToDestination[0][ lastNode ] > maxDepth ){
							break;
					
					}
				}

				if( reachAbleError && ( u != source ) && ( pathSize != 2 ) ) {
					continue;
				} 

				//everything is fine here 
				double total_nodescore = 0;
				double totalEdgeDist = 0;

				double prevNodeDist= scores[ i ][ k ].second;
				for(int pat = 1; pat < pathSize ; pat++ ) 
				{

					int pathNode = ex_edge[pat];
					double pres_edge = 0;
					for( int edgeIndex = 0; edgeIndex < edges[ ex_edge[ pat- 1 ] ].size(); edgeIndex++)
					{
						if( edges[ ex_edge[ pat - 1 ] ][ edgeIndex ] == ex_edge[ pat ] ) 
						{
							pres_edge = edgeWeight[ ex_edge[ pat - 1 ] ][ edgeIndex ];	
							break;
						}
					}

					totalEdgeDist += pres_edge;
					int offset= (DISCRETE_TIME * totalEdgeDist / delta_dist);
					double newDist = (prevNodeDist + totalEdgeDist);
					double timeOffset=0;
					if(TIME_VARYING)
					{
						timeOffset= ((int)(newDist/SPEED))/MINUTE_PER_HISTORY_SLOT;
					}

					double node_score = 0;

					if(USE_CLUSTER)
					{
						node_score = getClusterScore(pathNode, destinations, timeSlot + timeOffset, defaultAlpha, actualTravelDists, newDist, shortestPathDistances, 
						dataset, distanceToDestination, distanceFromDestination);
					}
					else
					{
						node_score = get_expected_trips(pathNode, destinations, timeSlot + timeOffset, defaultAlpha, actualTravelDists, newDist, shortestPathDistances, 
						dataset, distanceToDestination, distanceFromDestination);
					}
					// double node_score = getClusterScore(pathNode, destinations, timeSlot + ((int)(newDist/SPEED))/MINUTE_PER_HISTORY_SLOT, defaultAlpha, actualTravelDists, newDist, shortestPathDistances,
						// dataset, distanceToDestination, distanceFromDestination);
					total_nodescore += node_score;

				}

				// cout<<"huh3 "<<endl;

				double totaloffset1= (DISCRETE_TIME * totalEdgeDist / delta_dist);
				int newTime = (DISCRETE_TIME * (totalEdgeDist + prevNodeDist) / delta_dist);

				if(nodeToDagIndex[ lastNode ] > i && newTime <= DISCRETE_TIME)
				{
									
					if((scores[ nodeToDagIndex[ lastNode ] ][ newTime ].first  <  scores[ i ][ k ].first +  total_nodescore) || 
							((scores[ nodeToDagIndex[ lastNode ] ][ newTime ].first  ==  scores[ i ][ k ].first +  total_nodescore) && 
								scores[ i ][ k ].second + totalEdgeDist < scores[ nodeToDagIndex[ lastNode ] ][ newTime ].second))
					{
						scores[ nodeToDagIndex[ lastNode ] ][ newTime ].first = scores[ i ][ k ].first + total_nodescore;
						scores[ nodeToDagIndex[ lastNode ] ][ newTime ].second = scores[ i ][ k ].second + totalEdgeDist;
						store[ nodeToDagIndex[ lastNode ] ][ newTime ] = make_pair(ex_edge, k);
					}
				}

					
			}
		}
		
	}

	// for(int i=0; i<dag.size(); i++)
 //    {
 //        // vector <double> fogRow;
 //        int flag=0;
 //        for(int j=0; j<=DISCRETE_TIME; j++)
 //        {
 //             if(scores[i][j].first>=0)
 //             {
 //             	cout<<"i= "<<i<<" j= "<<j<<" sc= ";
 //             	cout<<scores[i][j].first<<" ";
 //             	flag=1;
 //             }
 //        }
 //        // scores.push_back(fogRow);
 //        if(flag==1)
 //        	cout<<endl;
 //    }

	// backtracking -> needed to change just like updation
	double max_score=-1;
	// long long int maxindex=-1;
	int max_time=-1;
	for( int popp=0;popp<=DISCRETE_TIME;popp++) 
	{
		//printf(" -- Score = %d (%.4f < %.4f?)\n",it->first, it->second.first, distFromSourceToDestination*alpha);
		if( scores[endIndex][popp].first>max_score ) 
		{
			max_score= scores[endIndex][popp].first;
			// maxindex=it->second.second.first;
			max_time= popp;
		}
	}

	cout<<max_score<<endl;
	cout<<max_time<<endl;

	if(max_score==-1)
	{
		vector <double> dummy(N);
		path= dijkstra_lengths(N, source, firstDestination, dummy, edges, edgeWeight);
		short_outputs+=1;
		cout<<"shortestpaths= "<<short_outputs<<endl;
		return make_pair(-1,path);
	}
	

	int traceLocation = endIndex;
	int tracetime = max_time;

	while (traceLocation != -1) 
	{
		// printf(" At %d: Time = %d (score: %.4f)\n",traceLocation, tracetime, scores[ traceLocation ][tracetime].first); fflush(stdout);
		path.push_back(dagIndexToNode[ traceLocation] );
		int prevtime = store[ traceLocation ][ tracetime ].second;
		vector<int> prev_locs = store[ traceLocation ][tracetime].first;
		for(int loc = ((int)prev_locs.size())-2; loc>=1; loc--)
		{
			path.push_back(prev_locs[loc]);
		}
		traceLocation = prev_locs[0];

		if(traceLocation < 0)
		{
			break;
		}
		traceLocation=nodeToDagIndex[traceLocation];

		tracetime = prevtime;
	}
	reverse(path.begin(), path.end());

	cout<<"BEST SCORE: dex 0.3  "<<max_score<<endl;
	double pathlen=getPathDist(path, edges, edgeWeight);
	cout<<"Pathlen= "<<pathlen<<endl;
	cout<<"delta_dist= "<<delta_dist<<endl;
	// for(int i=0;i<path.size();i++ )
	// {
	// 	cout<<path[i]<<endl;
	// }
	double pathlen2=distanceFromSource[firstDestination];
	cout<<"disjkstra Pathlen= "<<pathlen2<<endl;

	if (DEBUG && dag.size() <= 40) {
		printf("s:%d e:%d\n",startIndex, endIndex);
		for(int i = 0; i < dag.size(); i++) {
			int u = dag[i].second;
			double node_score = get_expected_trips(u, destinations, timeSlot, defaultAlpha, actualTravelDists, 0, shortestPathDistances, 
						dataset, distanceToDestination, distanceFromDestination);
			printf("%d %.4f\n",i,node_score);

		}
		for (int i = 0; i<dag.size(); i++) {
			int u = dag[i].second;
			for(int j = 0; j < edges[u].size(); j++) {
				int v = edges[u][j];
				if( nodeToDagIndex[ v ] == -1 )
					continue;
				//some nodes are adjacent to themselves in the graph
				if(u == v)
					continue;
				printf("%d %d %.4f\n",nodeToDagIndex[u], nodeToDagIndex[v],edgeWeight[ u ][ j ]);
			}
		}
		printf("Path: ");
		for (int i=0; i<path.size(); i++) {
			printf(" %d",nodeToDagIndex[path[i]]);
		}
		printf("\n");
	}

	return make_pair(max_score, path);
}
















