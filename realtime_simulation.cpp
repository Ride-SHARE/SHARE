#include <bits/stdc++.h>
#include <ctime>
#include "carpooling.h"
#include <chrono>

using namespace std;

// #define unordered_map map
// #define unordered_set set

char * timeSourceName, * locationInputName, * edgeInputName;
char * txtName;
string location;
string routeUsed,assignUsed;

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

vector< unordered_map<int, vector<request> > > sourceTimeDestination; 
vector< unordered_map<int, vector<request> > > passengerRequest;
int RUNNING_TIME=0;

/* Simulation */
int DAY = 6;
// 2013 Jan (21 is holiday)
int NY_WEEKDAYS[] = {2,3,4,7,8,9,10,11,14,15,16,17,18,22,23,24,25,28,29,30,31};
// 2012 Jan (23,24 is new year holiday)
int SG_WEEKDAYS[] = {3,4,5,6,9,10,11,12,13,16,17,18,19,20,25,26,27,30,31};
// 2008 Feb
int BJ_WEEKDAYS[] = {2,3,4,5,6,7,8};
// 2008 May-Jun
int SF_WEEKDAYS[] = {19,20,21,22,23,26,27,28,29,30,2,3,4,5,6,9,10};
set<int> weekdays;
int cabsShared = 0;
int globalPickups = 0, globalRequests = 0, globalRejections = 0;  
long long int passengerWaitTime[MAX_CAB_CAPACITY+1][MAX_END_TIME];
double globalDistTravelled[MAX_CAB_CAPACITY+1];
int countFalse = 0;
int assignPickups = 0;
int maxWaitedSlot = MINUTE_BEFORE_REJECTION/DELTA_TIME;
int queueTimeSlot = 0, remainingQueueSize = 0;
int total_count=0;
double percentageDiff=0;

vector< passenger > passengerPickedList;
// At source v - a list of passenger waiting for cab
unordered_map<int, vector< passenger > > passengerQueue;

/* Time analysis */
long long cnt_bestpath = 0, cnt_bestpath2 = 0;
double total_time_bestpath = 0, total_time_bestpath2 = 0;
clock_t clockStartTime[5];
//for runnning time
double optTotalTime1 = 0;
double optTotalTime2 = 0;
double optTotalTime2_1 = 0;
double optTotalTime2_2 = 0;
double optTotalTime2_3 = 0;
double optTotalTime2_4 = 0;
double optTotalTime2_5 = 0;
double optTotalTime2_6 = 0;
double optTotalTime2_7 = 0;
double optTotalTime2_8 = 0;
double optTotalTime2_9 = 0;

double optTotalTime3 = 0;
double optTotalTime3_1 = 0;
double optTotalTime3_2 = 0;
double optTotalTime3_3 = 0;
double optTotalTime3_4 = 0;
double optTotalTime3_5 = 0;
double optTotalTime3_6 = 0;
double optTotalTime3_7 = 0;
double optTotalTime3_8 = 0;
double optTotalTime3_9 = 0;
double optTotalTime4 = 0;


double dexTotalTime1 = 0;
double dexTotalTime2 = 0;
double dexTotalTime2_1 = 0;
double dexTotalTime2_2 = 0;
double dexTotalTime2_3 = 0;
double dexTotalTime2_4 = 0;
double dexTotalTime2_5 = 0;
double dexTotalTime2_6 = 0;
double dexTotalTime2_7 = 0;
double dexTotalTime2_8 = 0;
double dexTotalTime2_9 = 0;
double dexTotalTime3 = 0;
double dexTotalTime3_1 = 0;
double dexTotalTime3_2 = 0;
double dexTotalTime3_3 = 0;
double dexTotalTime3_4 = 0;
double dexTotalTime3_5 = 0;
double dexTotalTime3_6 = 0;
double dexTotalTime3_7 = 0;
double dexTotalTime3_8 = 0;
double dexTotalTime3_9 = 0;
double dexTotalTime4 = 0;

long numTrips1=0;
long numTrips2=0;
long numTrips2_1=0;
long numTrips2_2=0;
long numTrips2_3=0;
long numTrips2_4=0;
long numTrips2_5=0;
long numTrips2_6=0;
long numTrips2_7=0;
long numTrips2_8=0;
long numTrips2_9=0;

long numTrips3=0;
long numTrips3_1=0;
long numTrips3_2=0;
long numTrips3_3=0;
long numTrips3_4=0;
long numTrips3_5=0;
long numTrips3_6=0;
long numTrips3_7=0;
long numTrips3_8=0;
long numTrips3_9=0;
long numTrips4=0;

void start_clock(int slot = 0) {
	clockStartTime[slot] = clock();
}

double get_runtime(int slot = 0) {
	return float( clock()- clockStartTime[slot] )/ CLOCKS_PER_SEC;
}
/* Time analysis ends */

bool isReachable(int source, int destination) {
	vector< double > distanceFromSource( nodeID.size() );
	vector<int> path;
	path = dijkstra_lengths(nodeID.size(), source, destination, distanceFromSource, edges, edgeWeight);
	if(distanceFromSource[destination] == MAX_DISTANCE || !distanceFromSource[ destination ] ) {
		printf("Not reachable -\tSource: %d\tDestination: %d\n", source, destination);
		return false;
	}
	return true;
}

void updateTime(double duration, double durationDex, double shortestDistance)
{
	if(shortestDistance < 1)
	{
		optTotalTime1+= duration;
		dexTotalTime1 += durationDex;
		numTrips1+=1;
	}
	if(shortestDistance < 2){
		optTotalTime2 += duration;
		dexTotalTime2 += durationDex;
		numTrips2+=1;

	}
	if(shortestDistance < 2.1){
		optTotalTime2_1 += duration;
		dexTotalTime2_1 += durationDex;
		numTrips2_1+=1;

	}
	if(shortestDistance < 2.2){
		optTotalTime2_2 += duration;
		dexTotalTime2_2 += durationDex;
		numTrips2_2+=1;

	}
	if(shortestDistance < 2.3){
		optTotalTime2_3 += duration;
		dexTotalTime2_3 += durationDex;
		numTrips2_3+=1;

	}
	if(shortestDistance < 2.4){
		optTotalTime2_4 += duration;
		dexTotalTime2_4 += durationDex;
		numTrips2_4+=1;

	}
	if(shortestDistance < 2.5){
		optTotalTime2_5 += duration;
		dexTotalTime2_5 += durationDex;
		numTrips2_5+=1;

	}
	if(shortestDistance < 2.6){
		optTotalTime2_6 += duration;
		dexTotalTime2_6 += durationDex;
		numTrips2_6+=1;

	}
	if(shortestDistance < 2.7){
		optTotalTime2_7 += duration;
		dexTotalTime2_7 += durationDex;
		numTrips2_7+=1;

	}
	if(shortestDistance < 2.8){
		optTotalTime2_8 += duration;
		dexTotalTime2_8 += durationDex;
		numTrips2_8+=1;

	}
	if(shortestDistance < 2.9){
		optTotalTime2_9 += duration;
		dexTotalTime2_9 += durationDex;
		numTrips2_9+=1;

	}
	if(shortestDistance < 3){
		optTotalTime3 += duration;
		dexTotalTime3 += durationDex;
		numTrips3+=1;
		double optAvg3 = optTotalTime3/(double)numTrips3;
		cout<<"optAvgTime3= "<< optAvg3 << endl;
		double dexAvg3 = dexTotalTime3/(double)numTrips3;
		cout<<"dexAvgTime3= "<< dexAvg3 << endl;
	}
	if(shortestDistance < 3.1){
		optTotalTime3_1 += duration;
		dexTotalTime3_1 += durationDex;
		numTrips3_1+=1;

	}
	if(shortestDistance < 3.2){
		optTotalTime3_2 += duration;
		dexTotalTime3_2 += durationDex;
		numTrips3_2+=1;

	}
	if(shortestDistance < 3.3){
		optTotalTime3_3 += duration;
		dexTotalTime3_3 += durationDex;
		numTrips3_3+=1;

	}
	if(shortestDistance < 3.4){
		optTotalTime3_4 += duration;
		dexTotalTime3_4 += durationDex;
		numTrips3_4+=1;

	}
	if(shortestDistance < 3.5){
		optTotalTime3_5 += duration;
		dexTotalTime3_5 += durationDex;
		numTrips3_5+=1;

	}
	if(shortestDistance < 3.6){
		optTotalTime3_6 += duration;
		dexTotalTime3_6 += durationDex;
		numTrips3_6+=1;

	}
	if(shortestDistance < 3.7){
		optTotalTime3_7 += duration;
		dexTotalTime3_7 += durationDex;
		numTrips3_7+=1;

	}
	if(shortestDistance < 3.8){
		optTotalTime3_8 += duration;
		dexTotalTime3_8 += durationDex;
		numTrips3_8+=1;

	}
	if(shortestDistance < 3.9){
		optTotalTime3_9 += duration;
		dexTotalTime3_9 += durationDex;
		numTrips3_9+=1;

	}
	if(shortestDistance < 4){
		optTotalTime4 += duration;
		dexTotalTime4 += durationDex;
		numTrips4+=1;
	}

	cout<<endl;
	cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;
	cout<<percentageDiff/total_count<<endl;
}

vector<int> getBestPath(int source, vector<int> destinations, int timeSlot, double defaultAlpha, vector<double> &actualTravelDists, vector<double> &shortestPathDistances) {
	printf("Trip start for %d passengers\n", (int)destinations.size());

	printf("Source: %d(%lld)\tDestination:",source,nodeID[source]);
	for (int i=0; i<destinations.size(); i++) {
		printf(" %d(%lld)",destinations[i],nodeID[destinations[i]]);
	}
	printf("\tTime slot: %d\n",timeSlot); fflush(stdout);

	start_clock();

	int firstDestination = destinations[0];

	vector< double > distanceFromSource( nodeID.size() );
	vector<vector< double > > distanceFromDestination( destinations.size() );
	vector<vector< double > > distanceToDestination( destinations.size() );

	vector<int> dagExPath;
	dagExPath = dijkstra_lengths(nodeID.size(), source, firstDestination, distanceFromSource, edges, edgeWeight);

	if(distanceFromSource[firstDestination] == MAX_DISTANCE || !distanceFromSource[ firstDestination ] ) {
		printf("Cannot reach destination?\n");
		vector<int> emptyPath;
		return emptyPath;
	}

	for (int i=0; i<destinations.size(); i++) {
		distanceToDestination[i].resize( nodeID.size() );
		distanceFromDestination[i].resize( nodeID.size() );
		dijkstra_lengths(nodeID.size(), destinations[i], source, distanceToDestination[i], edgesReverse, edgeWeightReverse);
		dijkstra_lengths(nodeID.size(), destinations[i], source, distanceFromDestination[i], edges, edgeWeight);
	}

	if (ROUTE == 1) {
		printf("  Start DAG\n"); fflush(stdout);
		pair<double, vector<int> > pp;
		auto startDex = std::chrono::high_resolution_clock::now(); 
		pp = findDAGPath( nodeID.size(), source, destinations, timeSlot/MINUTE_PER_HISTORY_SLOT, defaultAlpha, actualTravelDists, shortestPathDistances,
			sourceTimeDestination, distanceFromSource, distanceToDestination, distanceFromDestination, edges, edgeWeight);
		
		auto stopDex = std::chrono::high_resolution_clock::now();
		auto durationDex = std::chrono::duration_cast<std::chrono::microseconds>(stopDex - startDex);

		double dagScore= pp.first;
		dagExPath= pp.second;

		if( RUN_OPTIMAL && shortestPathDistances[0] < 2 && dagScore!= -1 )
		{
			cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;

			auto start = std::chrono::high_resolution_clock::now(); 

			pair<double, vector<int> > optimal = optimalScore( nodeID.size(), source, destinations, timeSlot/MINUTE_PER_HISTORY_SLOT, defaultAlpha, actualTravelDists, shortestPathDistances,
				sourceTimeDestination, edges, edgeWeight, distanceFromSource, distanceToDestination, distanceFromDestination );
			auto stop = std::chrono::high_resolution_clock::now(); 
			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start); 

			vector<int> optimalPath= optimal.second;
			double optScore= optimal.first;

			cout<<"edges siz= "<<edges.size()<<endl;
			for(int pa=0;pa<optimalPath.size();pa++)
			{
				cout<<optimalPath[pa]<<" ";
			}
			if(optScore > 0)
			{

				total_count++;
				percentageDiff += (optScore-dagScore)/optScore;
				updateTime((double)duration.count(), (double)durationDex.count(), shortestPathDistances[0]);
				cout<<endl;
				cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;
				cout<<percentageDiff/total_count<<endl;
			}
		}

	}

	if (ROUTE == 2) {
		printf("  Start DAG EX\n"); fflush(stdout);
		pair<double, vector<int> > pp;

		cout<<"actDist= "<<actualTravelDists[0]<<endl<<endl;
		auto startDex = std::chrono::high_resolution_clock::now(); 

		pp = findDAGExtendedPath( nodeID.size(), source, destinations, timeSlot/MINUTE_PER_HISTORY_SLOT, defaultAlpha, actualTravelDists, shortestPathDistances,
			maxDepth, sourceTimeDestination, distanceFromSource, distanceToDestination, distanceFromDestination, edges, edgeWeight);
		auto stopDex = std::chrono::high_resolution_clock::now();
		auto durationDex = std::chrono::duration_cast<std::chrono::microseconds>(stopDex - startDex); 

		dagExPath = pp.second;
		double dexScore= pp.first;

		if( RUN_OPTIMAL && shortestPathDistances[0] < 3 && dexScore!= -1 )
		{
			cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;
			auto start = std::chrono::high_resolution_clock::now(); 
			
			pair<double, vector<int> > optimal = optimalScore( nodeID.size(), source, destinations, timeSlot/MINUTE_PER_HISTORY_SLOT, defaultAlpha, actualTravelDists, shortestPathDistances,
				sourceTimeDestination, edges, edgeWeight, distanceFromSource, distanceToDestination, distanceFromDestination );
			
			auto stop = std::chrono::high_resolution_clock::now(); 
			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start); 

			vector<int> optimalPath = optimal.second;
			double optScore = optimal.first;

			//cout<<"edges siz= "<<edges.size()<<endl;
			for(int pa=0;pa<optimalPath.size();pa++)
			{
				cout<<optimalPath[pa]<<" ";
			}

			if(optScore > 0)
			{
				total_count++;
				percentageDiff += (optScore-dexScore)/optScore;
				updateTime((double)duration.count(), (double)durationDex.count(), shortestPathDistances[0]);
			}
		}
	}

	if(ROUTE == 3 )
	{
		printf("  Start Optimal\n"); fflush(stdout);
		pair<double, vector<int> > pp;


		pp = findDAGPath( nodeID.size(), source, destinations, timeSlot/MINUTE_PER_HISTORY_SLOT, defaultAlpha, actualTravelDists, shortestPathDistances,
			 sourceTimeDestination, distanceFromSource, distanceToDestination, distanceFromDestination, edges, edgeWeight);

		dagExPath= pp.second;

		for(int pa=0;pa<dagExPath.size();pa++)
		{
			cout<<dagExPath[pa]<<" ";
		}
		cout<<endl;

		if(shortestPathDistances[0] < 3)
		{
			cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;

			pair<double, vector<int> > optimal  = optimalScore( nodeID.size(), source, destinations, timeSlot/MINUTE_PER_HISTORY_SLOT, defaultAlpha, actualTravelDists, shortestPathDistances,
				sourceTimeDestination, edges, edgeWeight, distanceFromSource, distanceToDestination, distanceFromDestination );
			
			vector<int> optimalPath = optimal.second;
			double optScore =  optimal.first;
			cout<<"edges siz= "<<edges.size()<<endl;
			for(int pa=0;pa<optimalPath.size();pa++)
			{
				cout<<optimalPath[pa]<<" ";
			}
			cout<<endl;
			cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;
		}

		//exit(1);
	}

	cnt_bestpath++;
	total_time_bestpath += get_runtime();
	cout<<"Trip Found"<<endl; fflush(stdout);
	
	// cout<<"CHECK ALPHA VIOLATION AMOUNT: "<<checkAlphaPath( dagExPath, edges, edgeWeight, distanceFromSource[ firstDestination ])<<endl;

	return dagExPath;
}

class cab {

private:
	int id;
	// Represents the path and current location of cab
	double fractionLastEdge; // How much fraction the cab has traversed in the last edge?
	vector< int > path;
	int lastNodeIndex; // index into the path array
	// Passenger information
	vector<passenger> passengerOnCab; // On cab already
	vector<passenger> passengerAssigned; // Not picked up yet
	vector<event> eventQueue; // e.g. o_1 pickup, o_2 pickup, o_1 dropoff, o_2 dropoff
	int numPickedPassenger; // Num of picked passenger in this round
	bool hasDroppedOff;
	bool hasPickedUp;

	double totalDistTravelled[ MAX_CAB_CAPACITY + 1]; // Total distance travelled with k passengers on cab

public:
	void initialize(passenger firstPass, vector<int> assignPath) {
		// TODO: reduce the score of the assign path so other cabs won't choose the same path
		id = firstPass.id;
		path = assignPath;
		lastNodeIndex = 0;
		fractionLastEdge = 0;
		passengerOnCab.clear();

		event es(firstPass.id, firstPass.source, false);
		eventQueue.push_back(es);
		event ed(firstPass.id, firstPass.destination, true);
		eventQueue.push_back(ed);
		addPassengerAssigned(firstPass);

		for (int i = 0; i <= cabCapacity; i++)
			totalDistTravelled[ i ] = 0;
	}

	void addPassenger( int currTimeSlot, passenger pass ) {
		pass.pickedTimeSlot = currTimeSlot;
		pass.passengerNo = passengerOnCab.size();
		pass.shortestPathDist = query(pass.source, pass.destination);
		pass.actualTravelDist = 0;

		passengerOnCab.push_back( pass );
		passengerPickedList.push_back( pass );

		if (passengerOnCab.size() > 1) {
			globalPickups++;
			for( int i = 0; i < passengerOnCab.size(); i++) {
				passengerOnCab[ i ].iShared = true;
			}
		}
	}

	void addPassengerAssigned(passenger pass) {
		pass.shortestPathDist = query(pass.source, pass.destination);
		passengerAssigned.push_back(pass);
	}

	void updateShortestPathFromEventQueue() {

		// The cab journey is over
		if (eventQueue.empty()) {
			return;
		}

		// Path should be modified to go from path[currNodeIndex] -> midStop -> dest
		int currNode = path[ lastNodeIndex ];
		// Erase old path
		path.erase(path.begin() + lastNodeIndex + 1, path.end());

		vector< double > dummy( nodeID.size() );
		int dest = eventQueue[0].node;
		vector< int > midPath = dijkstra_lengths(nodeID.size(), currNode, dest, dummy, edges, edgeWeight);
		path.insert(path.end(), midPath.begin() + 1, midPath.end());

		/*
		for( int i = 0; i < eventQueue.size() ; i++) {
			int dest = eventQueue[i].node;
			vector< int > midPath = dijkstra_lengths(nodeID.size(), currNode, dest, dummy, edges, edgeWeight);
			path.insert(path.end(), midPath.begin() + 1, midPath.end());
			currNode = dest;
		}
		*/

	}

	// This function is the same as updateShortestPathFromEventQueue if dijkstra algorithm is used
	void rerouteFirstPathFromEventQueue(int currTimeSlot) {

		printf("Calling reroute with eventQueue = %d\n",(int)eventQueue.size());

		// The cab journey is over
		if (eventQueue.empty()) {
			return;
		}

		// We don't allow the cab to detour if it is assigned a passenger but not yet pick up
		if (!passengerAssigned.empty()) {
			updateShortestPathFromEventQueue();
			return;
		}
		// We don't allow the cab to detour if the cab is full
		if (passengerOnCab.size() == cabCapacity) {
			updateShortestPathFromEventQueue();
			return;
		}

		printf("Original event queue\n");
		for (int i=0; i<eventQueue.size(); i++) {
			printf(" %d: %d at %d%c\n",i,eventQueue[i].pid, eventQueue[i].node, eventQueue[i].isDropoff ? 'D' : 'P');
		}
		fflush(stdout);

		int currNode = path[ lastNodeIndex ];
		if (currNode == eventQueue[0].node) {
			printf("Oops! Event didn't get processed!\n"); fflush(stdout);
		}


		vector< double > alphaV;
		vector< int > destList;
		vector<double> actualTravelDists,shortestPathDistances;
		double firstLegDist = query(currNode, eventQueue[0].node);
		double remainDist = 0;
		double minFirstLegAlpha = 10;
		for (int i=0; i<eventQueue.size(); i++) {
			if (i > 0) remainDist = query(eventQueue[i-1].node, eventQueue[i].node);

			int j;
			for (j=0; j<passengerOnCab.size(); j++) {
				if (passengerOnCab[j].id == eventQueue[i].pid) break;
			}
			
			double maxAllowDistance = alpha * passengerOnCab[j].shortestPathDist;
			double remainAlpha = (maxAllowDistance - passengerOnCab[j].actualTravelDist) / (firstLegDist + remainDist);
			double firstLegAlpha = (maxAllowDistance - remainDist - passengerOnCab[j].actualTravelDist) / firstLegDist;
			printf(" Alpha %d: %d first %.4f (%.4f/%.4f) remain %.4f (%.4f/%.4f)\n",i,eventQueue[i].pid,firstLegAlpha,(maxAllowDistance - remainDist - passengerOnCab[j].actualTravelDist),firstLegDist, remainAlpha,(maxAllowDistance - passengerOnCab[j].actualTravelDist), (firstLegDist + remainDist));
			alphaV.push_back(remainAlpha);
			minFirstLegAlpha = min(minFirstLegAlpha, firstLegAlpha);
			destList.push_back(passengerOnCab[j].destination);
			actualTravelDists.push_back(passengerOnCab[j].actualTravelDist);
			shortestPathDistances.push_back(passengerOnCab[j].shortestPathDist);
		}

		printf(" min first leg alpha = %.4f\n",minFirstLegAlpha);
		if (minFirstLegAlpha < 1.0) {
			printf("Oops! Min first leg alpha is wrong\n"); fflush(stdout);
			minFirstLegAlpha = 1.0001;
		}

		vector<int> bestPath;
		bestPath = getBestPath( currNode, destList, currTimeSlot, alpha, actualTravelDists, shortestPathDistances);
		// bestPath = getBestPath2( currNode, make_pair(destList[0], destList[1]), currTimeSlot, alphaV, minFirstLegAlpha).first;

		// Erase old path
		path.erase(path.begin() + lastNodeIndex + 1, path.end());

		path.insert(path.end(), bestPath.begin() + 1, bestPath.end());

		/*
		vector< double > dummy( nodeID.size() );
		for( int i = 1; i < eventQueue.size() ; i++) {
			int dest = eventQueue[i].node;
			vector< int > midPath = dijkstra_lengths(nodeID.size(), currNode, dest, dummy, edges, edgeWeight);
			path.insert(path.end(), midPath.begin() + 1, midPath.end());
			currNode = dest;
		}
		*/

		printf("Rerouting done\n"); fflush(stdout);
	}

	void processEvent(int currTimeSlot, event e) {
		printf("  Event happens at %d\n",e.node);
		if (e.isDropoff) {
			int j;
			for (j=0; j<passengerOnCab.size(); j++) {
				if (passengerOnCab[j].id == e.pid) break;
			}
			// keep track of passenger who shared
			if( passengerOnCab[ j ].iShared ) cabsShared += 1;

			printf("  Dropoff passenger %d[%d] at %d\n",passengerOnCab[j].passengerNo, passengerOnCab[j].id, path[ lastNodeIndex ]); fflush(stdout);
			if (passengerOnCab[j].actualTravelDist > alpha*(1.02)* passengerOnCab[ j ].shortestPathDist) {
				vector< double > distanceFromNode( nodeID.size() );
				dijkstra_lengths(nodeID.size(), passengerOnCab[ j ].source, passengerOnCab[ j ].destination, distanceFromNode, edges, edgeWeight);
				double dijDistance = distanceFromNode[ passengerOnCab[ j ].destination ];
				printf("    ASSERT! Alpha constraint violated: %.4f/%.4f = %.4f [Dij distance = %.4f]\n", passengerOnCab[j].actualTravelDist, passengerOnCab[j].shortestPathDist, passengerOnCab[j].actualTravelDist/passengerOnCab[j].shortestPathDist, dijDistance);
				countFalse += 1;
			}
			else 
				printf("    Alpha constraint is fine here: %.4f/%.4f = %.4f\n", passengerOnCab[j].actualTravelDist, passengerOnCab[j].shortestPathDist, passengerOnCab[j].actualTravelDist/passengerOnCab[j].shortestPathDist);
			passengerOnCab.erase( passengerOnCab.begin() + j);
			hasDroppedOff = true;
		}
		else {
			int j;
			for (j=0; j<passengerAssigned.size(); j++) {
				if (passengerAssigned[j].id == e.pid) break;
			}
			printf("  Pickup passenger %d[%d] at %d\n", (int)passengerOnCab.size(), passengerAssigned[j].id, path[ lastNodeIndex ]); fflush(stdout);
			addPassenger( currTimeSlot, passengerAssigned[j] );
			passengerAssigned.erase( passengerAssigned.begin() + j );
			hasPickedUp = true;
		}
	}

	bool simulate(int currTimeSlot, bool canAdhocPickup) {
		numPickedPassenger = 0;

		printf("Simulating cab%d at %d (%d passengers) - Source: %d - Destination: %d\n",id,path[lastNodeIndex], (int)passengerOnCab.size(), path[0], path[path.size()-1]); fflush(stdout);

		hasDroppedOff = false;
		hasPickedUp = false;
		while (!eventQueue.empty()) {
			if (eventQueue[0].node != path[lastNodeIndex]) break;

			printf(" Before simulation event\n");
			event e = eventQueue[0];
			processEvent(currTimeSlot, e);
			eventQueue.erase( eventQueue.begin() );
		}
		if (hasDroppedOff || hasPickedUp) {
			rerouteFirstPathFromEventQueue(currTimeSlot);
		}

		double timeLeft = DELTA_TIME;
		while (timeLeft > 0 && lastNodeIndex != path.size()-1) {
			// printf("   cab at %lld [%d/%d] with %d passenger, time left %.4f\n",path[lastNodeIndex],lastNodeIndex, (int)path.size(), (int)passengerOnCab.size(),timeLeft); fflush(stdout);
			double timeToNextNode = 0.0, distToNextNode = 0.0;
			for (int i = 0; i < edges[ path[lastNodeIndex] ].size(); i++) {
				if (edges[ path[lastNodeIndex] ][ i ] == path[ lastNodeIndex+1 ]) {
					timeToNextNode = edgeTime[ path[lastNodeIndex] ][ i ];
					distToNextNode = edgeWeight[ path[lastNodeIndex] ][ i ];
					break;
				}
				if (i == edges[ path[lastNodeIndex] ].size()-1) {
					printf("!!! ALERT AT %d to %d\n",path[lastNodeIndex],path[lastNodeIndex+1]); fflush(stdout);
				}
			}
			// printf("  Time / Weight %lld %lld -> %.4f / %.4f\n",path[ lastNodeIndex ],path[ lastNodeIndex + 1 ],timeToNextNode, distToNextNode); fflush(stdout);
			assert(lastNodeIndex < (path.size()-1));

 			if ( timeToNextNode * ( 1 - fractionLastEdge ) >= timeLeft ) {
				// Move car on partial edge
				fractionLastEdge += timeLeft / timeToNextNode ;
				timeLeft = 0;
			}
			else {
				// Add actual travel dist for each passenger
				for (int i = 0; i < passengerOnCab.size(); i++) {
					passengerOnCab[ i ].actualTravelDist += distToNextNode;
				}
				totalDistTravelled[ passengerOnCab.size() ] += distToNextNode;

				// Move car on full edge
				timeLeft -= timeToNextNode * ( 1 - fractionLastEdge );
				fractionLastEdge = 0;
				lastNodeIndex += 1;


				// Iterate on events
				hasDroppedOff = false;
				hasPickedUp = false;
				while (!eventQueue.empty()) {
					if (eventQueue[0].node != path[lastNodeIndex]) break;

					event e = eventQueue[0];
					processEvent(currTimeSlot, e);
					eventQueue.erase( eventQueue.begin() );
				}
				if (hasDroppedOff || hasPickedUp) {
					rerouteFirstPathFromEventQueue(currTimeSlot);
				}
				
				// Check if we can pick up another passenger here (only if no passenger is assigned)
				if ( lastNodeIndex != path.size()-1 && passengerOnCab.size() < cabCapacity && passengerAssigned.empty() ) {

					while (passengerOnCab.size() < cabCapacity) {
						int v = path[ lastNodeIndex ];
						vector<passenger> &passList = passengerQueue[ v ];
						int passengerFound = -1;
						pair<vector< event >, double> candPassenger;
						double minPassengerDist = MAX_DISTANCE;
						vector< event > midPointList;

						for (int i = 0; i < passList.size(); i++) {
							candPassenger = midPointCalculation(v, passengerOnCab, passList[i]);
							if( candPassenger.first.size() > 0 && candPassenger.second < minPassengerDist ) {
								minPassengerDist = candPassenger.second;
								midPointList = candPassenger.first;
								passengerFound = i;
								printf("    passenger cand %d(pid: %d): total dist = %.4f\n",i,passList[i].id,candPassenger.second);
							}
						}

						if( passengerFound >= 0 ) {
							assignPickups += 1;
							printf("  %dth passenger (pid: %d) picked up at %d [%d other passengers available]\n",(int)passengerOnCab.size()+1, passList[passengerFound].id, v, (int)passList.size()-1); fflush(stdout);
							addPassenger( currTimeSlot, passList[passengerFound] );
							passList.erase(passList.begin() + passengerFound);
							numPickedPassenger++;
							eventQueue = midPointList;
							rerouteFirstPathFromEventQueue(currTimeSlot);
						}
						else {
							break;
						}
					}
				}
				
			}

		}

		// This should be equivalent to eventQueue.empty()
		if ( lastNodeIndex == path.size()-1 ) {
			printf("  Cab arrives at destination!\n"); fflush(stdout);
			if (!eventQueue.empty()) {
				printf("Oops event queue is not empty: %d\n",(int)eventQueue.size());
				for (int i=0; i<eventQueue.size(); i++) {
					printf(" %d: %d at %d%c\n",i,eventQueue[i].pid, eventQueue[i].node, eventQueue[i].isDropoff ? 'D' : 'P');
				}
			}
		}

		// Arrived at destination?
		return (lastNodeIndex == path.size()-1);
	}

	int getNumPassenger() {
		return passengerOnCab.size()+passengerAssigned.size();
	}

	int getNode() {
		return path[ lastNodeIndex ];
	}

	vector<passenger> getPassengerOnCab() {
		return passengerOnCab;
	}

	vector<passenger> getPassengerAssigned() {
		return passengerAssigned;
	}

	int getNumPickedPassenger() {
		return numPickedPassenger;
	}

	int getDestination() {
		return path[ path.size()-1 ];
	}

	double getTotalDistTravelled(int k) {
		return totalDistTravelled[k];
	}

	vector< event > getEventQueue() {
		return eventQueue;
	}

	void setEventQueue(vector< event > newEventQueue) {
		eventQueue = newEventQueue;
	}

};

void getTimeSourceDestination() {
	ifstream file;
	file.open( timeSourceName );
	string s;
	sourceTimeDestination.resize( nodeID.size() );
	passengerRequest.resize( MAX_END_TIME / DELTA_TIME );
	int trips = 0;
	while( getline(file, s) ) 
	{
		stringstream ss( s );
		int date, source, timeSlot, dest;
		double rev = 0.0;
		ss>>date>>source>>timeSlot>>dest;
		if (location.compare("NY") == 0 || location.compare("SG") == 0) {
			ss >> rev;
		}
		

		timeSlot /= DELTA_TIME;
		// printf("Hello %d %d %d %d\n",date,source,timeSlot,dest);
		if (date > DAY) continue;
		if (weekdays.find(date) == weekdays.end()) continue;
		if (timeSlot >= endTime / DELTA_TIME) continue;

		if (source < nodeID.size() && dest < nodeID.size() && timeSlot >= startTime / DELTA_TIME) {
			trips++;

			request req;
			req.id = trips; req.source = source; req.timeSlot = timeSlot; req.destination = dest; req.revenue = rev;
			// printf("trips %d: D%d T%d %d(%.4f,%.4f)->%d(%.4f,%.4f)\n",trips,date,timeSlot,source,nodesToLatLon[source].first,nodesToLatLon[source].second,dest,nodesToLatLon[dest].first,nodesToLatLon[dest].second);

			if( date < DAY )
				sourceTimeDestination[ source ][ timeSlot/MINUTE_PER_HISTORY_SLOT ].push_back( req );
			else {
				passengerRequest[ timeSlot ][ source ].push_back( req );	
			}
		}
		
	}

	file.close();
	cout<<trips<<endl;
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
		
		index++;
		if(location.compare("BJ") == 0 && nodesToLatLon.size() == numNodes)
			break;
	}

	if (location.compare("SF") == 0 || location.compare("NY") == 0 || location.compare("SG") == 0) {
		Location.close();
		Location.open( edgeInputName );
	}

	cout<<"# of nodes= "<<nodeID.size()<<endl;
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

// Do On San Francisco

void takePartialGraphInput() {
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

		//set range
		double LatLimit = -122.435;
		double LonLimit = 37.7796;
		if(lat < LatLimit && lat > LatLimit - 0.01 && lon < LonLimit && lon > LonLimit - 0.01)
		{
			nodesToLatLon.push_back( make_pair( lat, lon) );
			nodeID.push_back( id );
			idToNode[ id ] = index;
			index++;
			//cout<<"aajao"<<endl;
		}

		if(location.compare("BJ") == 0 && nodesToLatLon.size() == numNodes)
			break;
	}

	if (location.compare("SF") == 0 || location.compare("NY") == 0 || location.compare("SG") == 0) {
		Location.close();
		Location.open( edgeInputName );
	}

	cout<<"Graph node number= "<<nodeID.size()<<endl<<endl;

	// Get edges
	int count = 0;
	edges.resize(nodeID.size());
	edgeWeight.resize(nodeID.size());
	edgeTime.resize(nodeID.size());
	edgesReverse.resize(nodeID.size());
	edgeWeightReverse.resize(nodeID.size());
	while( getline(Location, s) ) 
	{
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
		if( !oneWay ) 
		{
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


vector<int> initializeFreeCabs() {
	srand( 0 );
	vector<int> freeCabs;
	for (int i=0; i<availableCab; i++) {
		freeCabs.push_back(rand()%nodeID.size());
	}
	return freeCabs;
}

void addingPassenger(int timeSlot) {
	printf("== Adding passenger ==\n"); fflush(stdout);
	// Add all passengers who arrive at this timeSlot to passengerQueue
	for (unordered_map<int, vector<request> >::iterator iter = passengerRequest[ timeSlot ].begin(); 
				iter != passengerRequest[ timeSlot ].end(); iter++) {
		int u = iter->first;
		vector<request> requestList = iter->second;
		for (int i=0; i < requestList.size(); i++) {
			int v = requestList[ i ].destination;
			passenger p;
			p.id = requestList[i].id;
			p.startTimeSlot = timeSlot;
			p.pickedTimeSlot = -1;
			p.source = u;
			p.destination = v;
			p.iShared = false;
			p.revenue = requestList[i].revenue;
			if (isReachable(u,v)) {
				passengerQueue[ u ].push_back( p );
				remainingQueueSize++;
			}
		}			
	}
	printf("Remaining queue size = %d\n",remainingQueueSize);
}

struct candidate {
	int u,j;
	double value;
	vector<event> eq;
};

// Return the best delta distance given the current event queue and the new passenger
pair<double, vector<event> > getBestDeltaDistance( int timeSlot, cab c, passenger pass ) {
	int nodeOfCab = c.getNode();
	vector<event> eq = c.getEventQueue();
	vector<passenger> cp = c.getPassengerOnCab();
	vector<passenger> ap = c.getPassengerAssigned();

	int u = pass.source;
	int v = pass.destination;

	double shortestPathUV = query(u, v);
	
	event es(pass.id, u, false);
	event ed(pass.id, v, true);

	double bestDist = 1e10;
	vector<event> besteq;

	for (int sid = 0; sid < eq.size(); sid++) {
		for (int eid = sid; eid <= eq.size(); eid++) {
			// We will try to insert es into eq[sid] and ed into eq[eid], then check if this is feasible
			vector<event> candeq;
			candeq = eq;
			candeq.insert(candeq.begin()+eid, ed);
			candeq.insert(candeq.begin()+sid, es);

			unordered_map<int, double> totDist;
			unordered_map<int, double> shortDist;
			for (int k = 0; k < cp.size(); k++) {
				totDist[cp[k].id] = cp[k].actualTravelDist;
				shortDist[cp[k].id] = cp[k].shortestPathDist;
			}

			// Checking feasibility given candeq
			double sumDist = 0;
			bool feasible = true;
			for (int k = 0; k < candeq.size(); k++) {
				double segDist;
				if (k == 0) segDist = query(nodeOfCab, candeq[0].node);
				else segDist = query(candeq[k-1].node, candeq[k].node);
				sumDist += segDist;
				for (unordered_map<int, double>::iterator it=totDist.begin(); it!=totDist.end(); it++) {
					it->second += segDist;
				}
				if (candeq[k].isDropoff) {
					if (totDist[candeq[k].pid] > alpha * shortDist[candeq[k].pid]) {
						feasible = false; break;
					}
				}
				else {
					int aid;
					for (aid = 0; aid < ap.size(); aid++) {
						if (ap[aid].id == candeq[k].pid) break;
					}
					totDist[candeq[k].pid] = 0;
					if (aid < ap.size()) {
						shortDist[candeq[k].pid] = ap[aid].shortestPathDist;
						if (sumDist/SPEED + (timeSlot - ap[aid].startTimeSlot)*DELTA_TIME > MINUTE_BEFORE_REJECTION) {
							feasible = false; break;
						}
					}
					else {
						shortDist[candeq[k].pid] = shortestPathUV;
						if (sumDist/SPEED + (timeSlot - pass.startTimeSlot)*DELTA_TIME > MINUTE_BEFORE_REJECTION) {
							feasible = false; break;
						}
					}

					
				}
			}

			if (feasible && sumDist < bestDist) {
				bestDist = sumDist;
				besteq = candeq;
			}

		}
	}

	return make_pair(bestDist, besteq);
}

void requestPreliminaryCheck(int timeSlot, unordered_map<int, set< int > > &markDelete) {
	for (unordered_map<int, vector<passenger> >::iterator iter = passengerQueue.begin(); iter != passengerQueue.end(); iter++) {
		int u = iter->first;
		vector<passenger> passList = iter->second;

		for (int i = 0; i < passList.size(); i++) {
			// Rejection case
			if( (timeSlot - passList[ i ].startTimeSlot)*DELTA_TIME > MINUTE_BEFORE_REJECTION ) {
				markDelete[ u ].insert( i );
				globalRejections++;
				continue;
			}

			int v = passList[ i ].destination;

			// Not reachable case
			if ( !isReachable(u, v) ) {
				markDelete[ u ].insert( i );
				continue;
			}
		}
	}
}

void deleteProcessedRequests(unordered_map<int, set< int > > &markDelete, vector< int > &freeCabs, unordered_set<int> &deleteFreeCab) {
	printf("= Delete processed passenger requests =\n"); fflush(stdout);
	// Delete processed passenger requests
	for (unordered_map<int, set<int> >::iterator iter = markDelete.begin(); iter != markDelete.end(); iter++) {
		int u = iter->first;
		vector<passenger> passList = passengerQueue[ u ];
		set<int> toDelete = iter->second;
		vector<passenger> newPassList;

		for (int i=0; i < passList.size(); i++) {
			if (toDelete.find(i) == toDelete.end()) {
				newPassList.push_back(passList[i]);
			}
			else {
				remainingQueueSize--;
			}
		}
		passengerQueue[u] = newPassList;
	}

	printf("= Delete free cabs =\n"); fflush(stdout);
	vector<int> newFreeCabs;
	for (int i=0; i<freeCabs.size(); i++) {
		if (deleteFreeCab.find(i) == deleteFreeCab.end()) {
			newFreeCabs.push_back(freeCabs[i]);
		}
	}
	freeCabs = newFreeCabs;

}

void greedyPriceAssignment(int timeSlot, vector< cab > &cabs, vector< int > &freeCabs) {
	printf("== Start greedy price assignment ==\n"); fflush(stdout);
	start_clock(1);

	unordered_map<int, set< int > > markDelete;
	unordered_set<int> deleteFreeCab;

	requestPreliminaryCheck(timeSlot, markDelete);

	printf("Generating occupied cab pairs [Count = %d] (Time = %.4f)\n",(int)cabs.size(),get_runtime(1)); fflush(stdout);

	vector<vector<candidate> > pool(cabs.size()+freeCabs.size());
	// map pool index to <[free/occupy], id>
	vector<pair<bool, int> > poolmap(cabs.size()+freeCabs.size());

	for (int i = 0; i < cabs.size(); i++) {
		// Full cab, cannot assign more passenger
		if (cabs[i].getNumPassenger() == cabCapacity) {
			poolmap[i] = make_pair(false, i);
			continue;
		}

		int nodeOfCab = cabs[ i ].getNode();
		vector<event> eq = cabs[i].getEventQueue();

		double cabOldDist = 0.0;
		for (int k = 0; k < eq.size(); k++) {
			double segDist;
			if (k == 0) segDist = query(nodeOfCab, eq[0].node);
			else segDist = query(eq[k-1].node, eq[k].node);
			cabOldDist += segDist;
		}

		vector<candidate> vpool;

		for (unordered_map<int, vector<passenger> >::iterator iter = passengerQueue.begin(); iter != passengerQueue.end(); iter++) {
			int u = iter->first;
			vector<passenger> passList = iter->second;
			if (passList.empty()) continue;

			double distToPickup = query(nodeOfCab, u);
		
			for (int j = 0; j < passList.size(); j++) {
				// Impossible to pick up on time
				if (distToPickup/SPEED + (timeSlot - passList[j].startTimeSlot)*DELTA_TIME > MINUTE_BEFORE_REJECTION) continue;

				pair<double, vector<event> > bestPermutation = getBestDeltaDistance(timeSlot, cabs[i], passList[j]);
				double bestDist = bestPermutation.first;
				vector<event> besteq = bestPermutation.second;

				if (!besteq.empty()) {
					// (cabs[i], passList[j]) is a feasible order-vehicle pair
					double profit = passList[j].revenue - FUEL*(bestDist - cabOldDist);
					if (profit >= 0) {
						candidate pc;
						pc.u = u; pc.j = j; pc.value = profit; pc.eq = besteq;
						vpool.push_back(pc);
					}
				}
			}
		}

		pool[i] = vpool;
		poolmap[i] = make_pair(false, i);
	}

	printf("Generating free cab pairs [Count = %d] (Time = %.4f)\n",(int)freeCabs.size(),get_runtime(1)); fflush(stdout);

	for (int i = 0; i < freeCabs.size(); i++) {

		int nodeOfCab = freeCabs[ i ];

		vector<candidate> vpool;
		
		for (unordered_map<int, vector<passenger> >::iterator iter = passengerQueue.begin(); iter != passengerQueue.end(); iter++) {
			int u = iter->first;
			vector<passenger> passList = iter->second;
			if (passList.empty()) continue;

			double distToPickup = query(nodeOfCab, u);
		
			for (int j = 0; j < passList.size(); j++) {
				// Impossible to pick up on time
				if (distToPickup/SPEED + (timeSlot - passList[j].startTimeSlot)*DELTA_TIME > MINUTE_BEFORE_REJECTION) continue;

				int v = passList[j].destination;

				event es(passList[j].id, u, false);
				event ed(passList[j].id, v, true);

				vector<event> besteq;
				besteq.push_back(es);
				besteq.push_back(ed);

				double bestDist = distToPickup + query(u, v);

				// (cabs[i], passList[j]) is a feasible order-vehicle pair
				double profit = passList[j].revenue - FUEL*bestDist - BCOST;
				//printf("  cab %d - [%d,%d] profit = %.4f [%.4f - %.4f*%.4f - %.4f]\n",i,u,j,profit, passList[j].revenue, FUEL, bestDist, BCOST); fflush(stdout);
				if (profit >= 0) {
					candidate pc;
					pc.u = u; pc.j = j; pc.value = profit; pc.eq = besteq;
					vpool.push_back(pc);
				}
			}
		}

		pool[cabs.size()+i] = vpool;
		poolmap[cabs.size()+i] = make_pair(true, i);
	}

	printf("Generating done (Time = %.4f)\n",get_runtime(1)); fflush(stdout);

	double mxProfit = -1;
	pair<int,int> bestCab;

	do {
		/*
		printf("Pool size = %d\n",(int)pool.size());
		for (int i=0; i<pool.size(); i++) {
			printf("cab %c%d profit =",poolmap[i].first ? 'F' : 'O', poolmap[i].second);
			for (int j=0; j<pool[i].size(); j++) {
				printf(" %.2f", pool[i][j].profit);
			}
			printf("\n"); fflush(stdout);
		}
		*/

		mxProfit = -1;
		int countCandidate = 0;
		for (int i=0; i<pool.size(); i++) {
			for (int j=0; j<pool[i].size(); j++) {
				if (markDelete.find(pool[i][j].u) != markDelete.end() && markDelete[pool[i][j].u].find(pool[i][j].j) != markDelete[pool[i][j].u].end() )
					continue;
				countCandidate++;
				if (pool[i][j].value > mxProfit) {
					mxProfit = pool[i][j].value;
					bestCab = make_pair(i,j);
				}
			}
		}

		if (mxProfit < 0) break;

		int poolId = bestCab.first, candidateId = bestCab.second;
		candidate pc = pool[poolId][candidateId];

		printf(" Chosen [%d]! Profit = %.4f (%c%d)\n",countCandidate, mxProfit, poolmap[poolId].first ? 'F' : 'O', poolmap[poolId].second); fflush(stdout);

		bool isFreeCab = poolmap[poolId].first;
		if (isFreeCab) {
			int cabId = poolmap[poolId].second;

			// Initialize a cab
			cab freshcab;
			vector<int> path;
			path.push_back(freeCabs[poolmap[poolId].second]);
			freshcab.initialize(passengerQueue[pc.u][pc.j], path);
			freshcab.rerouteFirstPathFromEventQueue(timeSlot);

			globalRequests += 1;
			availableCab--;
			printf("Available Cab = %d\n",availableCab); fflush(stdout);
			deleteFreeCab.insert(poolmap[poolId].second);
			poolmap[poolId] = make_pair(false, cabs.size());
			cabs.push_back( freshcab );
		}
		else {
			int cabId = poolmap[poolId].second;
			// Carpool assign
			cabs[ cabId ].setEventQueue( pc.eq );
			cabs[ cabId ].addPassengerAssigned(passengerQueue[pc.u][pc.j]);
			cabs[ cabId ].rerouteFirstPathFromEventQueue(timeSlot);
		}

		// Update order
		markDelete[pc.u].insert(pc.j);

		// Update vehicle (now it must be having at least 1 passenger)
		int cabId = poolmap[poolId].second;
		int nodeOfCab = cabs[cabId].getNode();
		vector<event> eq = cabs[cabId].getEventQueue();

		// Full cab!
		if (cabs[cabId].getNumPassenger() == cabCapacity) {
			pool[poolId].clear();
			continue;
		}

		double cabOldDist = 0.0;
		for (int k = 0; k < eq.size(); k++) {
			double segDist;
			if (k == 0) segDist = query(nodeOfCab, eq[0].node);
			else segDist = query(eq[k-1].node, eq[k].node);
			cabOldDist += segDist;
		}

		for (int i=0; i<pool[poolId].size(); i++) {
			if (i == candidateId) continue; // Delete the order
			
			pool[poolId][i].value = -1; // Reset the profit of this (cab, vehicle) pair

			int nodeOfCab = cabs[cabId].getNode();
			vector<event> eq = cabs[cabId].getEventQueue();
			
			pair<double, vector<event> > bestPermutation = getBestDeltaDistance(timeSlot, cabs[cabId], passengerQueue[pc.u][pc.j]);
			double bestDist = bestPermutation.first;
			vector<event> besteq = bestPermutation.second;

			if (!besteq.empty()) {
				// (cabs[i], passList[j]) is a feasible order-vehicle pair
				double profit = passengerQueue[pc.u][pc.j].revenue - FUEL*(bestDist - cabOldDist);
				if (profit >= 0) {
					pc.value = profit; pc.eq = besteq;
					pool[poolId][i] = pc;
				}
			}
		}
	} while (mxProfit >= 0);

	printf("End of greedy price assignment (Time = %.4f)\n",get_runtime(1));

	deleteProcessedRequests(markDelete, freeCabs, deleteFreeCab);

}

void greedyShortestDeltaDistAssignment(int timeSlot, vector< cab > &cabs, vector< int > &freeCabs) {
	printf("== Start greedy shortest delta distance assignment ==\n"); fflush(stdout);
	start_clock(1);

	unordered_map<int, set< int > > markDelete;
	unordered_set<int> deleteFreeCab;

	requestPreliminaryCheck(timeSlot, markDelete);

	printf("Generating occupied cab pairs [Count = %d] (Time = %.4f)\n",(int)cabs.size(),get_runtime(1)); fflush(stdout);

	vector<vector<candidate> > pool(cabs.size()+freeCabs.size());
	// map pool index to <[free/occupy], id>
	vector<pair<bool, int> > poolmap(cabs.size()+freeCabs.size());

	for (int i = 0; i < cabs.size(); i++) {
		// Full cab, cannot assign more passenger
		if (cabs[i].getNumPassenger() == cabCapacity) {
			poolmap[i] = make_pair(false, i);
			continue;
		}

		int nodeOfCab = cabs[ i ].getNode();
		vector<event> eq = cabs[i].getEventQueue();

		double cabOldDist = 0.0;
		for (int k = 0; k < eq.size(); k++) {
			double segDist;
			if (k == 0) segDist = query(nodeOfCab, eq[0].node);
			else segDist = query(eq[k-1].node, eq[k].node);
			cabOldDist += segDist;
		}

		vector<candidate> vpool;

		for (unordered_map<int, vector<passenger> >::iterator iter = passengerQueue.begin(); iter != passengerQueue.end(); iter++) {
			int u = iter->first;
			vector<passenger> passList = iter->second;
			if (passList.empty()) continue;

			double distToPickup = query(nodeOfCab, u);
		
			for (int j = 0; j < passList.size(); j++) {
				// Impossible to pick up on time
				if (distToPickup/SPEED + (timeSlot - passList[j].startTimeSlot)*DELTA_TIME > MINUTE_BEFORE_REJECTION) continue;

				pair<double, vector<event> > bestPermutation = getBestDeltaDistance(timeSlot, cabs[i], passList[j]);
				double bestDist = bestPermutation.first;
				vector<event> besteq = bestPermutation.second;

				if (!besteq.empty()) {
					// (cabs[i], passList[j]) is a feasible order-vehicle pair
					double deltaDist = bestDist - cabOldDist;
					candidate pc;
					pc.u = u; pc.j = j; pc.value = deltaDist; pc.eq = besteq;
					vpool.push_back(pc);
				}
			}
		}

		pool[i] = vpool;
		poolmap[i] = make_pair(false, i);
	}

	printf("Generating free cab pairs [Count = %d] (Time = %.4f)\n",(int)freeCabs.size(),get_runtime(1)); fflush(stdout);

	for (int i = 0; i < freeCabs.size(); i++) {

		int nodeOfCab = freeCabs[ i ];

		vector<candidate> vpool;
		
		for (unordered_map<int, vector<passenger> >::iterator iter = passengerQueue.begin(); iter != passengerQueue.end(); iter++) {
			int u = iter->first;
			vector<passenger> passList = iter->second;
			if (passList.empty()) continue;

			double distToPickup = query(nodeOfCab, u);
		
			for (int j = 0; j < passList.size(); j++) {
				// Impossible to pick up on time
				if (distToPickup/SPEED + (timeSlot - passList[j].startTimeSlot)*DELTA_TIME > MINUTE_BEFORE_REJECTION) continue;

				int v = passList[j].destination;

				event es(passList[j].id, u, false);
				event ed(passList[j].id, v, true);

				vector<event> besteq;
				besteq.push_back(es);
				besteq.push_back(ed);

				double bestDist = distToPickup + query(u, v);

				// (cabs[i], passList[j]) is a feasible order-vehicle pair
				double deltaDist = bestDist;
				candidate pc;
				pc.u = u; pc.j = j; pc.value = deltaDist; pc.eq = besteq;
				vpool.push_back(pc);
			}
		}

		pool[cabs.size()+i] = vpool;
		poolmap[cabs.size()+i] = make_pair(true, i);
	}

	printf("Generating done (Time = %.4f)\n",get_runtime(1)); fflush(stdout);

	double minDeltaDist = 1e10;
	pair<int,int> bestCab;

	do {

		minDeltaDist = 1e10;
		int countCandidate = 0;
		for (int i=0; i<pool.size(); i++) {
			for (int j=0; j<pool[i].size(); j++) {
				if (markDelete.find(pool[i][j].u) != markDelete.end() && markDelete[pool[i][j].u].find(pool[i][j].j) != markDelete[pool[i][j].u].end() )
					continue;
				countCandidate++;
				if (pool[i][j].value < minDeltaDist ) {
					minDeltaDist = pool[i][j].value;
					bestCab = make_pair(i,j);
				}
			}
		}

		if (minDeltaDist >= 1e9) break;

		int poolId = bestCab.first, candidateId = bestCab.second;
		candidate pc = pool[poolId][candidateId];

		printf(" Chosen [%d]! Distance = %.4f (%c%d) [%d]\n",countCandidate, minDeltaDist, poolmap[poolId].first ? 'F' : 'O', poolmap[poolId].second, passengerQueue[pc.u][pc.j].id); fflush(stdout);

		bool isFreeCab = poolmap[poolId].first;
		if (isFreeCab) {
			int cabId = poolmap[poolId].second;

			// Initialize a cab
			cab freshcab;
			vector<int> path;
			path.push_back(freeCabs[poolmap[poolId].second]);
			freshcab.initialize(passengerQueue[pc.u][pc.j], path);
			freshcab.rerouteFirstPathFromEventQueue(timeSlot);

			globalRequests += 1;
			availableCab--;
			printf("Available Cab = %d\n",availableCab); fflush(stdout);
			deleteFreeCab.insert(poolmap[poolId].second);
			poolmap[poolId] = make_pair(false, cabs.size());
			cabs.push_back( freshcab );
		}
		else {
			int cabId = poolmap[poolId].second;
			// Carpool assign
			cabs[ cabId ].setEventQueue( pc.eq );
			cabs[ cabId ].addPassengerAssigned(passengerQueue[pc.u][pc.j]);
			cabs[ cabId ].rerouteFirstPathFromEventQueue(timeSlot);
		}

		// Update order
		markDelete[pc.u].insert(pc.j);

		// Update vehicle (now it must be having at least 1 passenger)
		int cabId = poolmap[poolId].second;
		int nodeOfCab = cabs[cabId].getNode();
		vector<event> eq = cabs[cabId].getEventQueue();

		// Full cab!
		if (cabs[cabId].getNumPassenger() == cabCapacity) {
			pool[poolId].clear();
			continue;
		}

		double cabOldDist = 0.0;
		for (int k = 0; k < eq.size(); k++) {
			double segDist;
			if (k == 0) segDist = query(nodeOfCab, eq[0].node);
			else segDist = query(eq[k-1].node, eq[k].node);
			cabOldDist += segDist;
		}

		for (int i=0; i<pool[poolId].size(); i++) {
			if (i == candidateId) continue; // Delete the order
			
			pool[poolId][i].value = 1e10; // Reset the profit of this (cab, vehicle) pair

			int nodeOfCab = cabs[cabId].getNode();
			vector<event> eq = cabs[cabId].getEventQueue();
			
			pair<double, vector<event> > bestPermutation = getBestDeltaDistance(timeSlot, cabs[cabId], passengerQueue[pc.u][pc.j]);
			double bestDist = bestPermutation.first;
			vector<event> besteq = bestPermutation.second;

			if (!besteq.empty()) {
				// (cabs[i], passList[j]) is a feasible order-vehicle pair
				double deltaDist = bestDist - cabOldDist;
				pc.value = deltaDist; pc.eq = besteq;
				pool[poolId][i] = pc;
			}
		}
	} while (minDeltaDist < 1e9);

	printf("End of greedy shortest delta distance assignment (Time = %.4f)\n",get_runtime(1));

	deleteProcessedRequests(markDelete, freeCabs, deleteFreeCab);

}

void greedyShortestPickupAssignment(int timeSlot, vector< cab > &cabs, vector< int > &freeCabs) {
	printf("== Start shortest pickup assignment ==\n"); fflush(stdout);
	start_clock(1);

	unordered_map<int, set< int > > markDelete;
	unordered_set<int> deleteFreeCab;

	requestPreliminaryCheck(timeSlot, markDelete);

	printf("Generating occupied cab pairs [Count = %d] (Time = %.4f)\n",(int)cabs.size(),get_runtime(1)); fflush(stdout);

	vector<vector<candidate> > pool(cabs.size()+freeCabs.size());
	// map pool index to <[free/occupy], id>
	vector<pair<bool, int> > poolmap(cabs.size()+freeCabs.size());

	for (int i = 0; i < cabs.size(); i++) {
		// Full cab, cannot assign more passenger
		if (cabs[i].getNumPassenger() == cabCapacity) {
			poolmap[i] = make_pair(false, i);
			continue;
		}

		int nodeOfCab = cabs[ i ].getNode();
		vector<event> eq = cabs[i].getEventQueue();

		vector<candidate> vpool;

		for (unordered_map<int, vector<passenger> >::iterator iter = passengerQueue.begin(); iter != passengerQueue.end(); iter++) {
			int u = iter->first;
			vector<passenger> passList = iter->second;
			if (passList.empty()) continue;

			double distToPickup = query(nodeOfCab, u);
			if (distToPickup/SPEED > MINUTE_BEFORE_REJECTION) continue;
		
			for (int j = 0; j < passList.size(); j++) {
				// Impossible to pick up on time
				if (distToPickup/SPEED + (timeSlot - passList[j].startTimeSlot)*DELTA_TIME > 5) continue;

				pair<double, vector<event> > bestPermutation = getBestDeltaDistance(timeSlot, cabs[i], passList[j]);
				vector<event> besteq = bestPermutation.second;
				double bestDist = 0;
				for (int k=0; k<besteq.size(); k++) {
					// Find out the distance for passList[j] pickup
					double segDist;
					if (k == 0) segDist = query(nodeOfCab, besteq[0].node);
					else segDist = query(besteq[k-1].node, besteq[k].node);
					bestDist += segDist;
					if (besteq[k].pid == passList[j].id) break;
				}

				if (!besteq.empty()) {
					// (cabs[i], passList[j]) is a feasible order-vehicle pair
					candidate pc;
					pc.u = u; pc.j = j; pc.value = bestDist; pc.eq = besteq;
					vpool.push_back(pc);
				}
			}
		}

		pool[i] = vpool;
		poolmap[i] = make_pair(false, i);
	}

	printf("Generating free cab pairs [Count = %d] (Time = %.4f)\n",(int)freeCabs.size(),get_runtime(1)); fflush(stdout);

	for (int i = 0; i < freeCabs.size(); i++) {

		int nodeOfCab = freeCabs[ i ];

		vector<candidate> vpool;
		
		for (unordered_map<int, vector<passenger> >::iterator iter = passengerQueue.begin(); iter != passengerQueue.end(); iter++) {
			int u = iter->first;
			vector<passenger> passList = iter->second;
			if (passList.empty()) continue;

			double distToPickup = query(nodeOfCab, u);
			if (distToPickup/SPEED > MINUTE_BEFORE_REJECTION) continue;
		
			for (int j = 0; j < passList.size(); j++) {
				// Impossible to pick up on time
				if (distToPickup/SPEED + (timeSlot - passList[j].startTimeSlot)*DELTA_TIME > MINUTE_BEFORE_REJECTION) continue;

				int v = passList[j].destination;

				event es(passList[j].id, u, false);
				event ed(passList[j].id, v, true);

				vector<event> besteq;
				besteq.push_back(es);
				besteq.push_back(ed);

				// (cabs[i], passList[j]) is a feasible order-vehicle pair
				candidate pc;
				pc.u = u; pc.j = j; pc.value = distToPickup; pc.eq = besteq;
				vpool.push_back(pc);
			}
		}

		pool[cabs.size()+i] = vpool;
		poolmap[cabs.size()+i] = make_pair(true, i);
	}

	printf("Generating done (Time = %.4f)\n",get_runtime(1)); fflush(stdout);

	double minDist = 1e10;
	pair<int,int> bestCab;

	do {

		minDist = 1e10;
		int countCandidate = 0;
		for (int i=0; i<pool.size(); i++) {
			for (int j=0; j<pool[i].size(); j++) {
				if (markDelete.find(pool[i][j].u) != markDelete.end() && markDelete[pool[i][j].u].find(pool[i][j].j) != markDelete[pool[i][j].u].end() )
					continue;
				countCandidate++;
				if (pool[i][j].value < minDist ) {
					minDist = pool[i][j].value;
					bestCab = make_pair(i,j);
				}
			}
		}

		if (minDist >= 1e9) break;

		int poolId = bestCab.first, candidateId = bestCab.second;
		candidate pc = pool[poolId][candidateId];

		printf(" Chosen [%d]! Distance = %.4f (%c%d) [%d]\n",countCandidate, minDist, poolmap[poolId].first ? 'F' : 'O', poolmap[poolId].second, passengerQueue[pc.u][pc.j].id); fflush(stdout);

		bool isFreeCab = poolmap[poolId].first;
		if (isFreeCab) {
			int cabId = poolmap[poolId].second;

			// Initialize a cab
			cab freshcab;
			vector<int> path;
			path.push_back(freeCabs[poolmap[poolId].second]);
			freshcab.initialize(passengerQueue[pc.u][pc.j], path);
			freshcab.rerouteFirstPathFromEventQueue(timeSlot);

			globalRequests += 1;
			availableCab--;
			printf("Available Cab = %d\n",availableCab); fflush(stdout);
			deleteFreeCab.insert(poolmap[poolId].second);
			poolmap[poolId] = make_pair(false, cabs.size());
			cabs.push_back( freshcab );
		}
		else {
			int cabId = poolmap[poolId].second;
			// Carpool assign
			cabs[ cabId ].setEventQueue( pc.eq );
			cabs[ cabId ].addPassengerAssigned(passengerQueue[pc.u][pc.j]);
			cabs[ cabId ].rerouteFirstPathFromEventQueue(timeSlot);
		}

		// Update order
		markDelete[pc.u].insert(pc.j);

		// Update vehicle (now it must be having at least 1 passenger)
		int cabId = poolmap[poolId].second;
		int nodeOfCab = cabs[cabId].getNode();
		vector<event> eq = cabs[cabId].getEventQueue();

		// Full cab!
		if (cabs[cabId].getNumPassenger() == cabCapacity) {
			pool[poolId].clear();
			continue;
		}

		for (int i=0; i<pool[poolId].size(); i++) {
			if (i == candidateId) continue; // Delete the order
			
			pool[poolId][i].value = 1e10; // Reset the value of this (cab, vehicle) pair

			int nodeOfCab = cabs[cabId].getNode();
			vector<event> eq = cabs[cabId].getEventQueue();
			
			pair<double, vector<event> > bestPermutation = getBestDeltaDistance(timeSlot, cabs[cabId], passengerQueue[pc.u][pc.j]);
			vector<event> besteq = bestPermutation.second;
			double bestDist = 0;
			for (int k=0; k<besteq.size(); k++) {
				// Find out the distance for passList[j] pickup
				double segDist;
				if (k == 0) segDist = query(nodeOfCab, besteq[0].node);
				else segDist = query(besteq[k-1].node, besteq[k].node);
				bestDist += segDist;
				if (besteq[k].pid == passengerQueue[pc.u][pc.j].id) break;
			}

			if (!besteq.empty()) {
				// (cabs[i], passList[j]) is a feasible order-vehicle pair
				pc.value = bestDist; pc.eq = besteq;
				pool[poolId][i] = pc;
			}
		}
	} while (minDist < 1e9);

	printf("End of greedy shortest pickup assignment (Time = %.4f)\n",get_runtime(1));

	deleteProcessedRequests(markDelete, freeCabs, deleteFreeCab);

	printf("== End shortest pickup assignment == (Time = %.4f)\n",get_runtime(1));

}

void matchingShortestEmptyAssignment(int timeSlot, vector< cab > &cabs, vector< int > &freeCabs, bool isFirstTimeSlot) {
	printf("== Matching shortest empty assignment ==\n"); fflush(stdout);
	while ( (availableCab > 0 ) && ( ( queueTimeSlot + MINUTE_BEFORE_ASSIGN_EMPTY_CAB/DELTA_TIME <= timeSlot ) || ( isFirstTimeSlot ) ) && ( queueTimeSlot <= timeSlot ) ) {

		unordered_map<int, set< int > > markDelete;
		vector< passenger > pickups;
		vector< pair<passenger, int> > pickupPassenger;

		printf("= queueTimeSlot : %d =\n", queueTimeSlot); fflush(stdout);

		for (unordered_map<int, vector<passenger> >::iterator iter = passengerQueue.begin(); 
					iter != passengerQueue.end(); iter++) {
			int u = iter->first;
			vector<passenger> passList = iter->second;
	
			for (int i = 0; i < passList.size(); i++) {
			
				if( passList[ i ].startTimeSlot > queueTimeSlot ) break;

				// Rejection case
				if( passList[ i ].startTimeSlot <= timeSlot - MINUTE_BEFORE_REJECTION/DELTA_TIME ) {
					markDelete[ u ].insert( i );
					globalRejections++;
					continue;
				}

				int v = passList[ i ].destination;

				// Not reachable case
				if ( !isReachable(u, v) ) {
					markDelete[ u ].insert( i );
					continue;
				}
			
				printf("Assign ID = %d (pid: %d)\n",(int)pickups.size(),passList[ i ].id);
				pickups.push_back( passList[i] ) ;
				pickupPassenger.push_back(make_pair( passList[ i ], i ) );
			}
		}
		// Hungarian Algorithm
		printf("= Start Hungarian =\n"); fflush(stdout);
		vector< int > cabLocation = freeCabs; 
		vector< int > cabIndex = freeCabs;
		vector< passenger > vp;
		vector< vector< passenger > > passengersOnBoard( freeCabs.size(), vp);
		
		for( int i = 0 ; i < cabs.size(); i++) {
			continue;
			int passengerCabCount = cabs[ i ].getNumPassenger( );
			if( passengerCabCount < cabCapacity ) {
				int nodeOfCab = cabs[ i ].getNode();
				// cab not busy
				if( nodeOfCab != -1 ) { 
					vector< passenger > passengers = cabs[ i ].getPassengerOnCab();
					cabLocation.push_back( nodeOfCab ); 
					passengersOnBoard.push_back( passengers );
					cabIndex.push_back(i);
				}
			}
		}

		vector< int > cabToPickupAssignment = assignFreeCabs( cabLocation, passengersOnBoard, pickups, timeSlot );
		printf("= End Hungarian =\n"); fflush(stdout);

		// Move available cabs to pick up assigned passengers
		vector< int > newFreeCabs;
		for (int i = 0; i < cabToPickupAssignment.size(); i++) {
			if (cabToPickupAssignment[ i ] == -1) {
				if( !passengersOnBoard[ i ].size() )
					newFreeCabs.push_back( freeCabs[ i ] );
				continue;
			}

			passenger pass = pickupPassenger[ cabToPickupAssignment[ i ] ].first;
			int v = pass.source;
			int w = pass.destination;
			
			markDelete[ v ].insert( pickupPassenger[ cabToPickupAssignment[i] ].second ); 

			int curPNum = passengersOnBoard[ i ].size();

			printf("  Assigned cab %d: [%d pass] from %d to %d (pid: %d)\n", i, curPNum, cabLocation[i], v, pass.id); fflush(stdout);
			
			if( curPNum == 0 ){
				// Free cab case

				// Initialize a cab
				cab freshcab;
				vector<int> path;
				path.push_back(cabLocation[ i ]);
				freshcab.initialize(pass, path);
				freshcab.rerouteFirstPathFromEventQueue(timeSlot);
				globalRequests += 1;
				availableCab--;
				printf("Available Cab = %d\n",availableCab);
				cabs.push_back( freshcab );
			} 
			else {
				// Pick up case
				vector< event > midPointList = midPointCalculation(cabLocation[ i ], passengersOnBoard[ i ], pass).first;

				if( midPointList.size() == 0 ) {
					cout<<"WTF WTF WTF"<<endl;
				}
				cabs[ cabIndex[ i ] ].setEventQueue( midPointList );
				cabs[ cabIndex[ i ] ].addPassengerAssigned(pass);
				cabs[ cabIndex[ i ] ].rerouteFirstPathFromEventQueue(timeSlot);
			} 
		
		}
		freeCabs = newFreeCabs;
		cout<<"FREE CABS SIZE: "<<freeCabs.size()<<endl;
		printf("= Delete processed passenger requests =\n"); fflush(stdout);
		// Delete processed passenger requests
		for (unordered_map<int, set<int> >::iterator iter = markDelete.begin(); iter != markDelete.end(); iter++) {
			int u = iter->first;
			vector<passenger> passList = passengerQueue[ u ];
			set<int> toDelete = iter->second;
			vector<passenger> newPassList;

			for (int i=0; i < passList.size(); i++) {
				if (toDelete.find(i) == toDelete.end()) {
					newPassList.push_back(passList[i]);
				}
				else {
					remainingQueueSize--;
				}
			}
			passengerQueue[u] = newPassList;
		}

		if (availableCab > 0) {
			queueTimeSlot++;
		}
	}

}

void simulateEachCab(int timeSlot, vector< cab > &cabs, vector< int > &freeCabs, bool canAdhocPickup ) {
	printf("== Simulate each existing cab ==\n"); fflush(stdout);
	// Simulate each existing cab
	vector< cab > newCabs;
	for (int i = 0; i < cabs.size(); i++) {
		bool terminate = cabs[ i ].simulate( timeSlot, canAdhocPickup );

		if ( cabs[ i ].getNumPickedPassenger() > 0 ) {
			// We picked up some passengers in this turn
			remainingQueueSize -= cabs[ i ].getNumPickedPassenger();
		}
		if (terminate) {
			// The cab arrives at destination
			availableCab++;
			freeCabs.push_back( cabs[ i ].getDestination() );
			// Update global dist calculation
			for (int j = 0; j <= cabCapacity; j++) {
				globalDistTravelled[j] += cabs[ i ].getTotalDistTravelled(j);
			}
		}
		else {
			newCabs.push_back( cabs[ i ] );
		}
	}

	cabs = newCabs;
}

void calculateMetrics(int timeSlot) {
	for (int i = 0; i < passengerPickedList.size(); i++) {
		if (passengerPickedList[ i ].pickedTimeSlot == timeSlot) {
			int waitedSlot = passengerPickedList[ i ].pickedTimeSlot - passengerPickedList[i].startTimeSlot;
			passengerWaitTime[passengerPickedList[ i ].passengerNo][ waitedSlot ]++;
			maxWaitedSlot = max(maxWaitedSlot, waitedSlot);
		}
	}
}

void printRunningTimes()
{

	double optAvgTime1 = optTotalTime1/(double)numTrips1;
	cout<<"optAvgTime1= "<< optAvgTime1 << endl;
	double dexAvgTime1 = dexTotalTime1/(double)numTrips1;
	cout<<"dexAvgTime1= "<< dexAvgTime1 << endl;

	double optAvgTime2 = optTotalTime2/(double)numTrips2;
	cout<<"optAvgTime2= "<< optAvgTime2 << endl;
	double dexAvgTime2 = dexTotalTime2/(double)numTrips2;
	cout<<"dexAvgTime2= "<< dexAvgTime2 << endl;

	double optAvgTime2_1 = optTotalTime2_1/(double)numTrips2_1;
	cout<<"optAvgTime2_1= "<< optAvgTime2_1 << endl;
	double dexAvgTime2_1 = dexTotalTime2_1/(double)numTrips2_1;
	cout<<"dexAvgTime2_1= "<< dexAvgTime2_1 << endl;

	double optAvgTime2_2 = optTotalTime2_2/(double)numTrips2_2;
	cout<<"optAvgTime2_2= "<< optAvgTime2_2 << endl;
	double dexAvgTime2_2 = dexTotalTime2_2/(double)numTrips2_2;
	cout<<"dexAvgTime2_2= "<< dexAvgTime2_2 << endl;

	double optAvgTime2_3 = optTotalTime2_3/(double)numTrips2_3;
	cout<<"optAvgTime2_3= "<< optAvgTime2_3 << endl;
	double dexAvgTime2_3 = dexTotalTime2_3/(double)numTrips2_3;
	cout<<"dexAvgTime2_3= "<< dexAvgTime2_3 << endl;

	double optAvgTime2_4 = optTotalTime2_4/(double)numTrips2_4;
	cout<<"optAvgTime2_4= "<< optAvgTime2_4 << endl;
	double dexAvgTime2_4 = dexTotalTime2_4/(double)numTrips2_4;
	cout<<"dexAvgTime2_4= "<< dexAvgTime2_4 << endl;

	double optAvgTime2_5 = optTotalTime2_5/(double)numTrips2_5;
	cout<<"optAvgTime2_5= "<< optAvgTime2_5 << endl;
	double dexAvgTime2_5 = dexTotalTime2_5/(double)numTrips2_5;
	cout<<"dexAvgTime2_5= "<< dexAvgTime2_5 << endl;

	double optAvgTime2_6 = optTotalTime2_6/(double)numTrips2_6;
	cout<<"optAvgTime2_6 = "<< optAvgTime2_6 << endl;
	double dexAvgTime2_6 = dexTotalTime2_6/(double)numTrips2_6;
	cout<<"dexAvgTime2_6 = "<< dexAvgTime2_6 << endl;

	double optAvgTime2_7 = optTotalTime2_7/(double)numTrips2_7;
	cout<<"optAvgTime2_7= "<< optAvgTime2_7 << endl;
	double dexAvgTime2_7 = dexTotalTime2_7/(double)numTrips2_7;
	cout<<"dexAvgTime2_7= "<< dexAvgTime2_7 << endl;

	double optAvgTime2_8 = optTotalTime2_8/(double)numTrips2_8;
	cout<<"optAvgTime2_8= "<< optAvgTime2_8 << endl;
	double dexAvgTime2_8 = dexTotalTime2_8/(double)numTrips2_8;
	cout<<"dexAvgTime2_8= "<< dexAvgTime2_8 << endl;

	double optAvgTime2_9 = optTotalTime2_9/(double)numTrips2_9;
	cout<<"optAvgTime2_9= "<< optAvgTime2_9 << endl;
	double dexAvgTime2_9 = dexTotalTime2_9/(double)numTrips2_9;
	cout<<"dexAvgTime2_9= "<< dexAvgTime2_9 << endl;
	

	double optAvgTime3 = optTotalTime3/(double)numTrips3;
	cout<<"optAvgTime3= "<< optAvgTime3 << endl;
	double dexAvgTime3 = dexTotalTime3/(double)numTrips3;
	cout<<"dexAvgTime3= "<< dexAvgTime3 << endl;

	double optAvgTime3_1 = optTotalTime3_1/(double)numTrips3_1;
	cout<<"optAvgTime3_1= "<< optAvgTime3_1 << endl;
	double dexAvgTime3_1 = dexTotalTime3_1/(double)numTrips3_1;
	cout<<"dexAvgTime3_1= "<< dexAvgTime3_1 << endl;

	double optAvgTime3_2 = optTotalTime3_2/(double)numTrips3_2;
	cout<<"optAvgTime3_2= "<< optAvgTime3_2 << endl;
	double dexAvgTime3_2 = dexTotalTime3_2/(double)numTrips3_2;
	cout<<"dexAvgTime3_2= "<< dexAvgTime3_2 << endl;

	double optAvgTime3_3 = optTotalTime3_3/(double)numTrips3_3;
	cout<<"optAvgTime3_3= "<< optAvgTime3_3 << endl;
	double dexAvgTime3_3 = dexTotalTime3_3/(double)numTrips3_3;
	cout<<"dexAvgTime3_3= "<< dexAvgTime3_3 << endl;

	double optAvgTime3_4 = optTotalTime3_4/(double)numTrips3_4;
	cout<<"optAvgTime3_4= "<< optAvgTime3_4 << endl;
	double dexAvgTime3_4 = dexTotalTime3_4/(double)numTrips3_4;
	cout<<"dexAvgTime3_4= "<< dexAvgTime3_4 << endl;

	double optAvgTime3_5 = optTotalTime3_5/(double)numTrips3_5;
	cout<<"optAvgTime3_5= "<< optAvgTime3_5 << endl;
	double dexAvgTime3_5 = dexTotalTime3_5/(double)numTrips3_5;
	cout<<"dexAvgTime3_5= "<< dexAvgTime3_5 << endl;

	double optAvgTime3_6 = optTotalTime3_6/(double)numTrips3_6;
	cout<<"optAvgTime3_6 = "<< optAvgTime3_6 << endl;
	double dexAvgTime3_6 = dexTotalTime3_6/(double)numTrips3_6;
	cout<<"dexAvgTime3_6 = "<< dexAvgTime3_6 << endl;

	double optAvgTime3_7 = optTotalTime3_7/(double)numTrips3_7;
	cout<<"optAvgTime3_7= "<< optAvgTime3_7 << endl;
	double dexAvgTime3_7 = dexTotalTime3_7/(double)numTrips3_7;
	cout<<"dexAvgTime3_7= "<< dexAvgTime3_7 << endl;

	double optAvgTime3_8 = optTotalTime3_8/(double)numTrips3_8;
	cout<<"optAvgTime3_8= "<< optAvgTime3_8 << endl;
	double dexAvgTime3_8 = dexTotalTime3_8/(double)numTrips3_8;
	cout<<"dexAvgTime3_8= "<< dexAvgTime3_8 << endl;

	double optAvgTime3_9 = optTotalTime3_9/(double)numTrips3_9;
	cout<<"optAvgTime3_9= "<< optAvgTime3_9 << endl;
	double dexAvgTime3_9 = dexTotalTime3_9/(double)numTrips3_9;
	cout<<"dexAvgTime3_9= "<< dexAvgTime3_9 << endl;
	

	double optAvgTime4 = optTotalTime4/(double)numTrips4;
	cout<<"optAvgTime4= "<< optAvgTime4 << endl;
	double dexAvgTime4 = dexTotalTime4/(double)numTrips4;
	cout<<"dexAvgTime4= "<< dexAvgTime4 << endl;

	return;
}

void printMetrics(int timeSlot, vector< cab > &cabs) {
	printf("Statistic for time slot (Stat %d)\n", timeSlot); fflush(stdout);
	printf("Total cabs dispatched: %d\n",globalRequests);
	printf("Total car pooling/shared cabs: %d %d\n",globalPickups, cabsShared);
	printf("Total rejection: %d\n",globalRejections);

	for (int k = 0; k <= maxWaitedSlot; k++) {
		printf("Waited %d =",k);
		for (int i = 0; i<cabCapacity; i++) {
			printf(" %lld",passengerWaitTime[ i ][ k ]);
		}
		printf("\n");
	}
	printf("Total distance travelled with k passengers\n");
	for (int k = 0; k <= cabCapacity; k++) {
		double baseDist = globalDistTravelled[ k ];
		for( int i = 0; i < cabs.size(); i++) {
			baseDist += cabs[ i ].getTotalDistTravelled( k );
		}
		printf("Passenger %d = %.4f\n", k, baseDist);
	}
	printf("TOTAL pickup on way: %d\n", assignPickups);
	printf("TOTAL ALPHA-Constraint violations: %d\n", countFalse);
	printf("Total best path called: %lld %lld\n", cnt_bestpath, cnt_bestpath2);
	printf("Total time used: %.4f %.4f\n", total_time_bestpath, total_time_bestpath2);
	printf("\n"); fflush(stdout);
}


int main(int argc, char const *argv	[])
{
	if (argc < 8) {
		printf("Usage: ./a.out [location={BJ|SF|NY|SG}] [alpha=1.3] [route={dij|dag|dex}] [maxDepth=0.2] [assign={hun|gsp|pxa|sdp}] [maxCab=2000] [cabCapacity=2] [Optional: starttime endtime]\n");
		return 0;
	}

	if (argc >= 10) {
		stringstream ssst(argv[8]);
		ssst>>startTime;
		stringstream sset(argv[9]);
		sset>>endTime;
	}

	assignUsed = argv[5];
	if (assignUsed.compare("hun") == 0) {
		ASSIGN = 0;
	}
	else if (assignUsed.compare("gsp") == 0) {
		ASSIGN = 1;
	}
	else if (assignUsed.compare("pxa") == 0) {
		ASSIGN = 2;
	}
	else if (assignUsed.compare("sdp") == 0) {
		ASSIGN = 3;
	}
	else {
		printf("Assignment algo not recognized.");
		return 1;
	}

	routeUsed = argv[3];
	if (routeUsed.compare("dij") == 0) {
		ROUTE = 0;
	}
	else if (routeUsed.compare("dag") == 0) {
		ROUTE = 1;
	}
	else if (routeUsed.compare("dex") == 0) {
		ROUTE = 2;
	}
	else if (routeUsed.compare("opt") == 0 )
	{
		ROUTE = 3;
	}
	else {
		printf("Route algo not recognized.");
		return 1;
	}

	printf("ASSIGN = %d\nROUTE = %d\n",ASSIGN,ROUTE); fflush(stdout);

	stringstream ssAva(argv[6]);
	ssAva>>availableCab;
	int maxCab = availableCab;

	stringstream ssCap(argv[7]);
	ssCap>>cabCapacity;
	printf("Available cabs = %d\nCab capacity = %d\n",availableCab,cabCapacity);

	location = argv[1];
	txtName = (char*)malloc(50);
	timeSourceName = (char*)malloc(50);
	locationInputName = (char*)malloc(50);
	edgeInputName = (char*)malloc(50);
	if (location.compare("BJ") == 0) {
		sprintf(txtName, "beijingIndex");
		sprintf(timeSourceName, "bj_output");
		sprintf(locationInputName, "bj_graph_time");
		weekdays.insert(BJ_WEEKDAYS, BJ_WEEKDAYS + sizeof(BJ_WEEKDAYS) / sizeof(int));
		DAY = 8;
	}
	else if (location.compare("NY") == 0) {
		sprintf(txtName, "nyIndex");
		sprintf(timeSourceName, "ny_output_price");
		sprintf(locationInputName, "ny_location");
		sprintf(edgeInputName, "ny_edge_time");
		weekdays.insert(NY_WEEKDAYS, NY_WEEKDAYS + sizeof(NY_WEEKDAYS) / sizeof(int));
		DAY = 24;
		BCOST = 2.5;
		FUEL = 1.43/60;
	}
	else if (location.compare("TOY") == 0) {
		sprintf(txtName, "toyIndex");
		sprintf(timeSourceName, "new_sf_output");
		sprintf(locationInputName, "new_sf_location");
		sprintf(edgeInputName, "new_sf_edge");
		weekdays.insert(SF_WEEKDAYS, SF_WEEKDAYS + sizeof(SF_WEEKDAYS) / sizeof(int));
		DAY = 28;
	}
	else if (location.compare("SF") == 0) {
		sprintf(txtName, "sfIndex");
		sprintf(timeSourceName, "sf_output");
		sprintf(locationInputName, "sf_location");
		sprintf(edgeInputName, "sf_edge");
		weekdays.insert(SF_WEEKDAYS, SF_WEEKDAYS + sizeof(SF_WEEKDAYS) / sizeof(int));
		DAY = 28;
	}
	else if (location.compare("SG") == 0) {
		sprintf(txtName, "sgIndex");
		sprintf(timeSourceName, "sg_output_price");
		sprintf(locationInputName, "sg_location");
		sprintf(edgeInputName, "sg_edge");
		weekdays.insert(SG_WEEKDAYS, SG_WEEKDAYS + sizeof(SG_WEEKDAYS) / sizeof(int));
		DAY = 19;
		BCOST = 2.8;
		FUEL = 0.22; // Based on the gas price (SGD) at 2012, assume 1 litre can travel 10km
	}
	else {
		printf("Location not recognized.");
		return 1;
	}
	
	stringstream ssParam(argv[4]);	
	stringstream ssAlpha(argv[2]);
	ssParam>>maxDepth;
	ssAlpha>>alpha;
	printf("Max Depth = %.2f\nAlpha = %.2f\n", maxDepth, alpha); fflush(stdout);

	takeGraphInput();
	printf("Finish taking graph input\n"); fflush(stdout);
	getTimeSourceDestination();
	printf("Finish get time source destination\n"); fflush(stdout);

	queryInit( txtName );
	intializeFrequent(nodeID.size(), FREQUENCY_OPT, sourceTimeDestination, edges, edgeWeight);
	cout<<"abbe"<<endl;
	assignExtendEdge(nodeID.size(), maxDepth, edges, edgeWeight);
	cout<<"okay"<<endl;
	generateNearbyNodes(nodeID.size(), edges, edgeWeight);

	// MAIN SIMULATION
	memset(passengerWaitTime, 0, sizeof(passengerWaitTime));
	queueTimeSlot = startTime/DELTA_TIME;
	bool canAdhocPickup = true;
	vector< cab > cabs;
	vector< int > freeCabs = initializeFreeCabs();
	for (int k = 0; k <= cabCapacity; k++) {
		globalDistTravelled[k] = 0;
	}

	for(int globalTime = startTime; globalTime < endTime || availableCab != maxCab; globalTime += DELTA_TIME) {
		
		int timeSlot = globalTime / DELTA_TIME;

		printf("=== Starting at time slot %d ===\n",timeSlot); fflush(stdout);
		if (globalTime >= MAX_END_TIME) {
			printf("Max end time reached %d -- FORCE STOP\n",MAX_END_TIME);
			break;
		}
		printf("Remaining queue size = %d\n",remainingQueueSize);

		addingPassenger(timeSlot);

		bool isFirstTimeSlot;
		switch (ASSIGN) {
			case 0 : 
				isFirstTimeSlot = (globalTime == startTime);
				matchingShortestEmptyAssignment(timeSlot, cabs, freeCabs, isFirstTimeSlot);
				break;
			case 1 :
				greedyShortestPickupAssignment(timeSlot, cabs, freeCabs);
				break;
			case 2 :
				greedyPriceAssignment(timeSlot, cabs, freeCabs);
				break;
			case 3 :
				greedyShortestDeltaDistAssignment(timeSlot, cabs, freeCabs);
				break;
		}

		simulateEachCab(timeSlot, cabs, freeCabs, canAdhocPickup);

		calculateMetrics(timeSlot);
		printMetrics(timeSlot, cabs);

		if(RUN_OPTIMAL)
			printRunningTimes();


	
	}

	
	

	if (RUN_OPTIMAL) {
		double netPercent= percentageDiff*100/total_count;
		cout<<"netPercent= "<< netPercent << endl;
	}

	return 0;
}





