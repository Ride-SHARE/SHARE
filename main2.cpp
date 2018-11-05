#include <bits/stdc++.h>
#include <ctime>
#include "ioD.h"

using namespace std;

long long n;
int queryCnt;
char * txtName, * degName, * outName, *inName;
char * timeSourceName, * locationInputName, * edgeInputName;
edgeL * deg;
edgeS * labelout, *labelin;
edgeS * labelx, * labely;
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

// Takes in node index
double query(int x, int y)
{
	// if we already have this key pair value then return 
	if( frequentPickup[ x ] && frequentDrop[ y ] ) {
		return distanceFrequentNodes[ x ][ y ];
	}
	
	if( timeOptimize.find( make_pair(x, y) ) != timeOptimize.end() ) {
		return timeOptimize[ make_pair(x, y) ]; 
	}

	if (x == y) return 0;
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
	timeOptimize[ make_pair(x, y) ] = float(ans * locationFactor)/1000;

	// For NY it is 2, otherwise it is 1
	return float(ans * locationFactor)/1000;
}

void loadIndex()
{
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
// to mark the nodes in the network
vector<int> marked;
int DEBUG = 1;
bool PRINT_PATH = false;
int TOTAL_TRIP = 30;
float alpha  = 1.3;
double maxDepth = 0.2;
long long int beta = 0;
long long int maxDistance = 100000;
/* nodes global structure */
vector< pair<double, double> > nodes;
vector<long long int> nodeID;
map<long long int, long long int> idToNode;
map< long long int, pair<double, double> > nodeToLatLon;

/* Edges global structure */
vector< vector<long long int> > edges;
vector< vector<double> > edgeWeight;
vector< vector<long long int> > edgesReverse;
vector< vector<double> > edgeWeightReverse;

vector< map<long long int, vector<long long int> > > sourceTimeDestination; 

int countWin = 0;

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

// TODO: To get rid of
double distFromSourceToDestination; 

void getTimeSourceDestination_Other() {
	ifstream file;
	file.open( timeSourceName );
	string s;
	int trips = 0;
	sourceTimeDestination.resize(nodeID.size());
	while( getline(file, s) ) {
		stringstream ss( s );
		long long int source, timeSlot, dest;
		ss>>source>>timeSlot;
		source = idToNode[source];
		while( ss>>dest ) {
			trips++;
			dest = idToNode[dest];
			sourceTimeDestination[ source ][ timeSlot ].push_back( dest );
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
	int trips = 0;
	while( getline(file, s) ) {
		stringstream ss( s );
		long long int source, timeSlot, dest;
		ss>>source>>timeSlot>>dest;
		if (source < nodeID.size() && dest < nodeID.size()) {
			trips++;
			sourceTimeDestination[ source ][ timeSlot ].push_back( dest );
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

vector<long long int>  dijkstra_lengths(long long int S, long long int D, vector< double > &distanceFromSource,
	vector< vector<long long int> > &edges, vector< vector<double> > &edgeWeight) { 
	
   	vector<long long int> prevNode(nodeID.size());
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

	while(dj.size()) { 
		x = dj.top();
		dj.pop();

		u = x.second;
		for(int i=0; i< edges[ u ].size(); i++) { 	
			v = edges[ u ][ i ];
			if( !marked[ v ])
				continue;
			alt = distanceFromSource[ u ] + edgeWeight[ u ][ i ];
			if(alt < distanceFromSource[ v ])
			{ 	distanceFromSource[ v ] = alt;
				dj.push( make_pair(-alt, v) );
				prevNode[v] = u;
			}
		}
	}


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

void printPath(vector<long long int> &path, vector< long long int > &weights) {
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
		printf("  %7lld:\t%lld\t(Dist: %.4f) [%.6f]\n",nodeID[v], weights[v], cumDist, nodeToLatLon[ nodeID[ v ] ].first);
	}
	fflush(stdout);
}

int getPathScore(vector<long long int> &path, vector< long long int > &weights) {
	int score = 0;
	for(int i=1;i<path.size()-1; i++) {
		score += weights[ path[i] ];
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

vector<long long int> HD(long long int S, long long int D, vector< long long int > &expectedTrips) { 
	vector<float> dist(nodeID.size());
 
	for(int i=0; i< nodeID.size() ; i++)
	{ 	dist[ i ] = 100000;
	}

	dist[S] = 0;

	// ( (score, node), (distance, path) )
	priority_queue< pair< pair<float, long long int>, pair< float, vector<long long int> > > > dj;
	
	vector<long long int> path;
	path.push_back(S);
	dj.push( make_pair( make_pair(0, S), make_pair(0, path) ) );

	while(dj.size())
	{ 
		pair< pair<float, long long int>, pair< float, vector<long long int> >  > x = dj.top();
		dj.pop();
		
		long long int u = x.first.second;
		
		if( x.second.first > distFromSourceToDestination )
			continue; 
		
		if( u == D ) {
			return x.second.second;
		}
		
		for(int i=0; i< edges[u].size(); i++)
		{ 	long long int v = edges[u][i];
			if( !marked[ v ] )
				continue;

			/* BUG: should skip the whole i loop 
			for(int j=0; j<x.second.second.size(); j++) 
				if( x.second.second[j] == v)
					continue;
			*/
			if (find(x.second.second.begin(), x.second.second.end(), v) != x.second.second.end()) {
				continue;
			}

			float alt = -x.first.first + (edgeWeight[u][i]/(expectedTrips[v]+1) );
			if( (x.second.first + edgeWeight[u][i]) < dist[v] )
			{ 	dist[v] = (x.second.first + edgeWeight[u][i]); 
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


void get_all_trips(long long int source,long long int start_time, vector<long long int> &tripList) {

	if( source < sourceTimeDestination.size() ) {
		if( sourceTimeDestination[source].find(start_time) != sourceTimeDestination[source].end() ) {

			for(int i = 0; i < sourceTimeDestination[ source][start_time].size(); i++) {
				tripList.push_back( sourceTimeDestination[source][ start_time ][ i ] );
			}

		}
	}
}

long long int get_expected_trips(long long int source, long long int destination, long long int midStop,
 	long long int start_time) {  	

	vector<long long int> tripList;
	
	get_all_trips(midStop, start_time, tripList);

  	int counter = 0;

  	long long int s = source;
  	long long int d = destination;
  	long long int v = midStop;
	double tvd, tsv;
	// no need ot calculate them for every iteration 
	if( tripList.size() ) {
		tsv = distanceFromSource[ v ];
  		tvd = distanceToDestination[ v ];
  	}

  	for(int i = 0; i < tripList.size(); i++) {

  		long long int w = tripList[ i ];
 		double tvw = query(v, w);
		double tdw = distanceFromDestination[ w ];

	  	if( tvd + tdw < alpha*tvw ) {
	  		counter += 1;
	  	}
	  	else {

	  		double twd = distanceToDestination[ w ];
	  		double distanceA = ( tsv + tvw + twd);
	  		double distanceB = alpha*distFromSourceToDestination;
	  		if( distanceA <= distanceB) {
	  			counter += 1;
	  		}
	  	}
  	}

  	return counter;
}

vector<long long int> dagScore1(long long int source, long long int destination, vector< double > &distanceFromSource,
	vector< double > &distanceToDestination, vector< long long int > &weights) {
	vector< pair< double, long long int> > dag;

	//cout<<distanceFromSource[ destination ]<<endl;
	for(int i=0; i <  nodes.size() ; i++) {	
		if( !marked[ i ] )
			continue;

		if( ( distanceFromSource[i] + distanceToDestination[i] ) < alpha*distFromSourceToDestination ) {
			dag.push_back( make_pair(-distanceToDestination[ i ], i ) );
		}
	}

	sort( dag.begin(), dag.end() );
	// (Distance, Last node) at dag[i] given a particular score
	vector< map<int, pair<double, long long int> > > scores( dag.size() );
  
  	vector< long long int> nodeToDagIndex(nodeID.size());
	// cout<<dag.size()<<endl;
	for (int i=0; i<nodeID.size(); i++) {
  		nodeToDagIndex[i] = -1;
  	}
	for(int i=0; i <  dag.size() ; i++) {
		nodeToDagIndex[ dag[i].second ] = i;
	}
		
	long long int startIndex = nodeToDagIndex[source];
  	long long int endIndex = nodeToDagIndex[destination];

  	scores[ startIndex ][ 0 ] = make_pair(0.0, -1);
  	for(int i = startIndex; i < dag.size(); i++) {
  		//printf("At index %d size %d [%lld -> %lld]\n",i,(int)scores[ i ].size(),startIndex,endIndex);
  		// Maintain the monocity of scores[i]
		double lastDist = 1e10; // infinity    	
		vector<int> delScore;
    	for(map<int, pair<double, long long int> >::reverse_iterator it=scores[ i ].rbegin(); it != scores[ i ].rend(); it++) { 
    		if (it->second.first >= lastDist) {
    			delScore.push_back(it->first);
    			continue;
    		}
    		lastDist = it->second.first;
    	}
    	for (int j = 0; j < delScore.size(); j++) {
    		scores[i].erase(delScore[j]);
    	}

    	long long int u = dag[i].second;
    	// cout<<weights[u]<<endl;
	    for(int j = 0; j < edges[u].size(); j++) {
	    	long long int v = edges[u][j];
	      	long long int vIndex = nodeToDagIndex[ v ];
	      	if( nodeToDagIndex[ v ] == -1 )
	        	continue;
	      	//some nodes are ajacent to themselves in the graph
	      	if(u == v)
	        	continue;

	      	if( vIndex >= i ) {
	      		for(map<int, pair<double, long long int> >::iterator it=scores[ i ].begin(); it != scores[ i ].end(); it++)
	        	{ 
	        		// cout<<it->first<<endl;
	        		if( it->first > 10000)
	        			continue;

	        		int curScore = it->first + weights[v];
	        		double prevDist = it->second.first;
	          		if( scores[ nodeToDagIndex[v] ].find( curScore ) != scores[ nodeToDagIndex[v] ].end() ) {
			          	if( scores[ nodeToDagIndex[v] ][curScore].first >= prevDist + edgeWeight[u][j] ) {
			            	scores[ nodeToDagIndex[v] ][curScore] = make_pair(prevDist + edgeWeight[u][j], u);
			            }
			        }
			        else
			        {	scores[ nodeToDagIndex[v] ][curScore] = make_pair(prevDist + edgeWeight[u][j], u);	
			        }
				}
			}
		}
	}
	long long int bestScore = 0;

	for( map<int, pair<double, long long int> >::iterator it=scores[ endIndex ].begin(); it != scores[ endIndex ].end(); it++) {
		//printf(" -- Score = %d (%.4f < %.4f?)\n",it->first, it->second.first, distFromSourceToDestination*alpha);
		if( it->second.first < distFromSourceToDestination*alpha ) {
			bestScore = it->first;
		}
	}
	vector<long long int> path;
	long long int traceLocation = destination;
	int traceScore = bestScore;
	while (traceLocation != -1) {
		//printf(" At %lld[%lld]: Score = %d (dist: %.4f)\n",traceLocation, nodeToDagIndex[traceLocation], traceScore, scores[ nodeToDagIndex[traceLocation] ][traceScore].first);
		path.push_back(traceLocation);
		int prevScore = traceScore - weights[traceLocation];
		traceLocation = scores[ nodeToDagIndex[traceLocation] ][traceScore].second;	
		traceScore = prevScore;
	}
	reverse(path.begin(), path.end());

	//cout<<"BEST SCORE:   "<<bestScore<<endl;
	return path;
}


void extendEdges(long long int source, long long int node, vector<long long int> &path, vector< pair<double, vector<long long int> > >  &paths,
	double pathLen) {
	
	for( int i = 0; i < path.size(); i++) {
		if( node == path[i] )
			return;
	}
	path.push_back( node );

	if( path.size() > 1 ) {
		paths.push_back( make_pair(pathLen, path) ) ;
	}

	if( pathLen >= maxDepth ) {
		path.pop_back();
		return ;
	}
	
	for(long long int j=0; j < edges[ node ].size(); j++) {
		long long int newNode = edges[ node ][j];
		if( !marked[ newNode ] )
			continue;
		extendEdges( source, newNode, path, paths, pathLen + edgeWeight[ node ][j] ) ;
	}

	path.pop_back();
}

vector<long long int> dagExtendedEdges(long long int source, long long int destination, vector< double > &distanceFromSource, vector< double > &distanceToDestination,
  vector< vector< pair<double, vector<long long int> > > > &extendEdge, vector< vector< long long int > > &extendEdgeWeights, vector< long long int > &weights) {

	vector< pair< double, int> > dag;
	//printf("Start extended %.2f\n",( clock() - startTimeTrip ) / (float) CLOCKS_PER_SEC);

	for(int i=0; i <  nodes.size() ; i++) {	
		if( !marked[ i ] )
			continue;
		if( ( distanceFromSource[i] + distanceToDestination[i] ) < alpha*distFromSourceToDestination ) {
			dag.push_back( make_pair(-distanceToDestination[ i ], i ) );
		}
	}

	sort( dag.begin(), dag.end() );

	// (Distance, (Last Node, Ex Edge used) ) at dag[i] given a particular score
  	vector< map<int, pair<double, pair<int,int> > >  > scores(dag.size());
  	
  	vector< int> nodeToDagIndex(nodeID.size());
  	for (int i=0; i<nodeID.size(); i++) {
  		nodeToDagIndex[i] = -1;
  	}
	for(int i=0; i <  dag.size() ; i++) {
		nodeToDagIndex[ dag[i].second ] = i;
	}
	
	int startIndex = nodeToDagIndex[ source ];
  	int endIndex = nodeToDagIndex[ destination ];

  	scores[ startIndex ][ 0 ] = make_pair(0, make_pair(-1,0)) ;
  	vector< set<int> > reach(nodeID.size());
  	vector< set<int> > reachFrom(nodeID.size());

  	//printf("Start iterate %.2f\n",( clock() - startTimeTrip ) / (float) CLOCKS_PER_SEC);

  	int countDead = 0;
  	int count_used = 0,count_reached = 0,count_edge = 0,count_relax = 0,count_update = 0,count_ops1 = 0,count_ops2 = 0,count_get_trip = 0;
  	for(int i = startIndex; i < dag.size(); i++) {
		int u = dag[i].second;
		
		// Maintain the monocity of scores[i]
		double lastDist = 1e10; // infinity		
		vector<int> delScore;
		for(map<int, pair<double, pair<int,int> > >::reverse_iterator it=scores[ i ].rbegin(); it != scores[ i ].rend(); it++) { 
			if (it->second.first >= lastDist) {
				delScore.push_back(it->first);
				continue;
			}
			lastDist = it->second.first;
		}
		for (int j = 0; j < delScore.size(); j++) {
			scores[i].erase(delScore[j]);
		}

		// We need not explore this DAG node
		if ( scores[i].empty() ) {
			countDead += 1;
			continue;
		}
		
		for(int j = 0; j < extendEdge[u].size(); j++) {
			
			pair<double, vector<long long int> > v = extendEdge[u][j];
			bool continueLoop = false;
			bool invalidNode = false;
			bool reachAbleError = false;
			int pathSize = v.second.size();

			count_edge++;

			if( pathSize < 1)
				continue;

			int lastNode =  v.second[ pathSize - 1 ];

			if ( lastNode == u )
				continue;

			// All nodes in the extended path except the last node has to be before u in the DAG
			for(int k = 1; k < pathSize ; k++ ) {
				int pathNode = v.second[ k]; 
				
				if( ( nodeToDagIndex[ pathNode ] <= i ) &&  ( k == ( pathSize -1 ) ) ) {
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

			// Avoid cycle while taking reverse path
			for(int k1 = 0; k1 < pathSize ; k1++ ) {	
				for(int k2 = k1 + 1; k2 < pathSize; k2++ ) {
					if( reach[ v.second[ k2 ] ].find( v.second[ k1 ] ) != reach[ v.second[ k2 ] ].end() ) {
						reachAbleError = true;
					}
				}
			}

			count_reached++;

		  	if( reachAbleError ) {
				continue;
		  	}

			count_used++;

		  	// Update reach and reachFrom so that we won't have cycles in the path
		  	for(int k1 = 0; k1 < pathSize ; k1++ ) {
				for(int k2 = k1 + 1; k2 < pathSize; k2++ ) {
					reach[ v.second[ k1 ] ].insert( v.second[ k2 ] );
					reachFrom[ v.second[ k2] ].insert( v.second[ k1 ] );
				}
			} 

			for(int k1 = pathSize - 1; k1 >= 0 ; k1-- ) {
				/*
				for( set<int>::iterator sit = reachFrom[ v.second[ k1 ] ].begin(); sit != reachFrom[ v.second[ k1 ] ].end(); sit++) {
					reach[ (*sit) ].insert( v.second.begin()+k1+1, v.second.end() );
					count_ops1++;
				}
				for(int k2 = k1 + 1; k2 < pathSize; k2++ ) {
					reachFrom[ v.second[ k2] ].insert( reachFrom[ v.second[ k1 ] ].begin(), reachFrom[ v.second[ k1 ] ].end() );
					count_ops2++;
				}
				*/
				for( set<int>::iterator sit = reachFrom[ v.second[ k1 ] ].begin(); sit != reachFrom[ v.second[ k1 ] ].end(); sit++) {
					for(int k2 = k1 + 1; k2 < pathSize; k2++ ) {
						// Two checks are needed because it is a directed graph
						if (abs(distanceFromSource[(*sit)] - distanceFromSource[v.second[k2]]) > maxDepth && abs(distanceToDestination[(*sit)] - distanceToDestination[v.second[k2]]) > maxDepth)
							continue;
						reach[ (*sit) ].insert( v.second[k2] );
						reachFrom[ v.second[ k2] ].insert( (*sit) );
					}
				}
				
			}

			// The extended edge weight is calculated here (lazy calculation)
			int vWeight = 0;
			for(int k = 1; k < pathSize ; k++ ) {
				int pathNode = v.second[ k];
				if (weights[pathNode] == -1) {
					count_get_trip++;
					weights[pathNode] = get_expected_trips(source, destination, pathNode, timeSlot, sourceTimeDestination);
				}
				vWeight += weights[pathNode];
			}

			int vIndex = nodeToDagIndex[ lastNode ];
		  	extendEdgeWeights[u][j] = vWeight;

	  		for(map<int, pair<double, pair<int,int> > >::iterator it=scores[ i ].begin(); it != scores[ i ].end(); it++) { 

	  			count_relax++;

			  	if( it->first > 10000)
						continue;

				int curScore = it->first + vWeight;
				double prevDist = it->second.first;
				double curDist = v.first;
		  		if( scores[ vIndex ].find( curScore ) != scores[ vIndex ].end() ) {
				  	if( scores[ vIndex][curScore].first >= prevDist + curDist ) {
						scores[ vIndex ][curScore] = make_pair(prevDist + curDist, make_pair(u,j));;
					}
				}
				else {
					scores[ vIndex ][curScore] = make_pair(prevDist + curDist, make_pair(u,j));;	
				}
			}
		}
		
	}
	cout<<"  DEAD EX :"<< countDead<<" "<<dag.size()- startIndex <<endl;
	long long int bestScore = 0;

	for( map<int, pair<double, pair<int,int> > >::iterator it=scores[ endIndex ].begin(); it != scores[ endIndex ].end(); it++) {	
		if( it->second.first < distFromSourceToDestination*alpha ) {
			bestScore = it->first;
		}
	}

	//printf("Start tracing %.2f\n",( clock() - startTimeTrip ) / (float) CLOCKS_PER_SEC);

	vector<long long int> path;
	int traceLocation = destination;
	int traceScore = bestScore;
	while (traceLocation != -1) {
		//if (traceScore > 0) printf(" At %d[%d]: Score = %d (dist: %.4f)\n",traceLocation, nodeToDagIndex[traceLocation], traceScore, scores[ nodeToDagIndex[traceLocation] ][traceScore].first); fflush(stdout);
		pair<int,int> prevInfo = scores[ nodeToDagIndex[traceLocation] ][traceScore].second;
		if (prevInfo.first != -1) {
			//if (extendEdgeWeights[prevInfo.first][prevInfo.second] > 0) {
			//	printf(" At %d[%lld]: Score = %d (this score = %lld)\n",traceLocation, nodeID[traceLocation], traceScore, extendEdgeWeights[prevInfo.first][prevInfo.second]); fflush(stdout);
			//}
			vector<long long int> extendPath = extendEdge[prevInfo.first][prevInfo.second].second;
			for (int i=extendPath.size()-1; i>0; i--) {
				path.push_back(extendPath[i]);
			}
		}
		else {
			path.push_back(traceLocation);
		}

		traceLocation = prevInfo.first;	
		if (traceLocation != -1) {
			traceScore -= extendEdgeWeights[prevInfo.first][prevInfo.second];
		}
	}
	reverse(path.begin(), path.end());
	//printf("Edge = %d\tReached = %d\tUsed = %d\tRelax = %d\tUpdate = %d (%d,%d)\n",count_edge,count_reached,count_used,count_relax,count_update,count_ops1,count_ops2);
	//printf("Get trip = %d\tTotal possible trip = %d\tTheo score = %lld\n",count_get_trip,count_valid_dag_trip,bestScore);
	//printf("Path size = %d\n",(int)path.size());

	return path;
}

long long int calculateScore( vector< long long int > &path, vector< long long int> &expectedTrips) {
	long long int ret = 0;
	for(int i = 1; i < path.size(); i++) {
		ret += expectedTrips[ path[i] ];
	}
	return ret;
}


int bfs(int source, int destination, map<int, int> &path, vector<long long int> &scores,
	int score, double pathLen, double pathLenMax) {
	
	if( ( source == destination ) && ( pathLen <= pathLenMax ) ) {
		return score - scores[ destination ];
	}
	
	if( ( clock() - startTimeTrip ) / (float) CLOCKS_PER_SEC > 60) {
		return -1;
	}

	if( ( pathLen + query(source, destination) ) >= pathLenMax + 0.05 )
		return 0;

	path[ source ] = 1;
	int maxx = -1;
	//cout<<( pathLen + query(source, destination) ) << " "<<pathLenMax<<endl;

	for(int j = 0; j < edges[ source ].size(); j++) {
		int newNode = edges[ source ][ j ];
		if( path[ newNode ] )
			continue;
		if( !marked[ newNode] )
				continue;
		maxx = max( maxx, bfs( newNode, destination, path, scores,
		score + scores[ newNode ], pathLen + edgeWeight[ source ][ j ] , pathLenMax) );
	}	

	path[ source ] = 0;

	return maxx;
}

bool potentialScore(long long int source, long long int destination, long long int timeSlot, 
	vector< vector< pair<double, vector<long long int> > > > &extendEdge)
{
	cout<<"Trip start"<<endl;
	printf("Source: %lld(%lld)\tDestination: %lld(%lld)\tTime slot: %lld\n",source,nodeID[source],destination,nodeID[destination],timeSlot);
	cout<<endl;
	distanceFromSource.resize( nodeID.size() );
	distanceFromDestination.resize( nodeID.size() );
	distanceToDestination.resize( nodeID.size() );

	vector<long long int> dijkstraPath = dijkstra_lengths(source, destination, distanceFromSource, edges, edgeWeight);
	dijkstraTime[ tripNumber ] = float( clock()- startTimeTrip )/ CLOCKS_PER_SEC;

	dijkstra_lengths(destination, source, distanceFromDestination, edges, edgeWeight);
	
	if(distanceFromSource[destination] == maxDistance || !distanceFromSource[ destination ] || ( query(source, destination ) > 4 ) ) {
		printf("Cannot reach destination?\n");
		return false;
	}
	
	dijkstra_lengths(destination, source, distanceToDestination, edgesReverse, edgeWeightReverse);

	if( DEBUG ) {
		printf("  [Dijkstra length] %.4f\n",distanceFromSource[ destination ]);
	}

	vector<long long int> expectedTrips( nodeID.size() ); 
	
	// From s to d, via v, where v is within 1+alpha
	for( int i =0; i<nodeID.size(); i++) {	
		if( !marked[ i ] )
			continue;
		if( ( distanceFromSource[ i ] + distanceToDestination[ i ] ) 
				< alpha	*distanceFromSource[ destination ] ) {
			//long long int passengersNumber;
			expectedTrips[ i ] = get_expected_trips(source, destination, i, timeSlot);
		}
	}

	// cache this value
	distFromSourceToDestination = query(source, destination );
	printf("  [Query length] %.4f\n",distFromSourceToDestination);

	/* Dijkstra path */
	startTimeTrip = clock();
	map<int, int> path;
	int dijkstraPathScore = bfs(source, destination, path, expectedTrips, 0, 0, distFromSourceToDestination*alpha );
	dijkstraTime[ tripNumber ] = ( clock()- startTimeTrip ) / (float) CLOCKS_PER_SEC;
	dijkstraScore[ tripNumber ] = dijkstraPathScore+1;
	dijkstraDist[ tripNumber ] = getPathDist(dijkstraPath);
	printf("  [Dijkstra score] %d (Time: %.2f)\n",dijkstraPathScore  + 1,dijkstraTime[ tripNumber ]);
	if (PRINT_PATH) {
		printPath(dijkstraPath, expectedTrips);
	}

	/* Max Score / Distance path */
	startTimeTrip = clock();
	distFromSourceToDestination = distanceFromSource[ destination ];
	
	vector<long long int> maxScorePerDistancePath = HD( source, destination, expectedTrips);
	int maxScorePerDistancePathScore = getPathScore(maxScorePerDistancePath, expectedTrips);
	maxScorePerDistanceTime[ tripNumber ] = ( clock() - startTimeTrip ) / (float) CLOCKS_PER_SEC;
	maxScorePerDistanceScore[ tripNumber ] = maxScorePerDistancePathScore+1;
	maxScorePerDistanceDist[ tripNumber ] = getPathDist(maxScorePerDistancePath);
	printf("  [Max Score/Distance path score] %d (Time: %.2f)\n",maxScorePerDistancePathScore +1,maxScorePerDistanceTime[ tripNumber ] );
	if (PRINT_PATH) {
		printPath(maxScorePerDistancePath, expectedTrips);
	}


	/* calculating DAG parameters */
	startTimeTrip = clock();
	vector<long long int> dagPath = dagScore1(source, destination, distanceFromSource,
	distanceToDestination, expectedTrips);
	int scoreDag = getPathScore(dagPath, expectedTrips);
	
	//cout<<scoreMinPath<<" SCORES "<<scoreDag<<endl;
	dagTimeTaken[ tripNumber ] = float( clock()- startTimeTrip )/ CLOCKS_PER_SEC;
	printf("  [DAG score] %d (Time: %.2f)\n",scoreDag+1,dagTimeTaken[ tripNumber ]);
	if (PRINT_PATH) {
		printPath(dagPath, expectedTrips);
	}

	dagScore[ tripNumber ] = scoreDag + 1 ;
	dagDist[ tripNumber ] = getPathDist(dagPath);

	/* calculating backDAG	 parameters */
	startTimeTrip = clock();
	
	vector< vector< long long int> > extendEdgeWeights(nodeID.size());
	for( int i = 0; i < nodeID.size(); i++) {		
		if( !marked[ i ] )
			continue;
		for(int j =0; j < extendEdge[i].size(); j++) {
			long long int scoreAlongPath = calculateScore( extendEdge[i][j].second, expectedTrips);
			extendEdgeWeights[ i ].push_back( scoreAlongPath );
		}
	}

	vector<long long int> dagExPath = dagExtendedEdges( source, destination, distanceFromSource, distanceToDestination, extendEdge, extendEdgeWeights, expectedTrips);
	int scoreExDag = getPathScore(dagExPath, expectedTrips);
	dagExScore[ tripNumber ] = scoreExDag + 1 ;
	dagExTime[ tripNumber ] = ( clock() - startTimeTrip ) / (float) CLOCKS_PER_SEC;
	dagExDist[ tripNumber ] = getPathDist(dagExPath);
	printf("  [Extended DAG score] %d (Time: %.2f)\n",scoreExDag+1,dagExTime[ tripNumber ]);
	if (PRINT_PATH) {
		printPath(dagExPath, expectedTrips);
	}
	timeOptimize.clear();
	//cout<<	"SCORE ALONG EX "<<scoreDag<<"  RATIO WITH DAG "<< float( scoreDag + 1)/(dagScore[ tripNumber ]*minPathScore[ tripNumber] )<<endl;
	return true;	
}


/* new function: we store distance to and from nodes which are in 
top 10 percentile wrt number of trips from them
*/
bool intializeFrequent() {
	// frequent pickups
	vector<int> v;
	for( int i = 0; i< nodeID.size(); i++ ) {
		int count = 0;
		for( int j = 0; j < 96; j++) {
			count += sourceTimeDestination[ i ][ j ].size() ; 
		}
		v.push_back( count );
	}
	sort( v.begin(), v.end() );
	int threshold = v[ v.size()*90/100 ];
	
	for( int i = 0; i< nodeID.size(); i++ ) {
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
	vector<int> des( nodeID.size(), 0);
	for( int i = 0; i< nodeID.size(); i++ ) {
		for( int j = 0; j < 96; j++) {
			for( int k = 0; k < sourceTimeDestination[ i ][ j ].size(); k++) {
				des[ sourceTimeDestination[ i ][ j ][ k ] ] += 1;
			} 
		}
	}

	vector< int > desSort = des; 
	sort( desSort.begin(), desSort.end() );
	threshold = desSort[ desSort.size()*90/100 ];
	
	for( int i = 0; i< nodeID.size(); i++ ) {
		if(  des[ i ] >= threshold ) {
			frequentDrop[ i ] = 1;
		}
		else 
			frequentDrop[ i ] = 0; 
	}

	for( int i = 0; i< nodeID.size(); i++ ) {	
		if( !frequentPickup[ i ] ) 
			continue;

	//	cout<<i<<" "<<nodeID.size()<<endl;
		vector< double > distanceFromSourceL( nodeID.size() );
		dijkstra_lengths(i, -1, distanceFromSourceL, edges, edgeWeight);
	
		for( int j = 0;  j < nodeID.size(); j++ ) {
			if( frequentDrop[ j ] ) {
				distanceFrequentNodes[ i ][ j ] = distanceFromSourceL[ j ]; 
			} 
		}
		//cout<<distanceFrequentNodes[i].size()<<endl;
	}
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

	degName = (char*)malloc(1+strlen(txtName) + 50);
	sprintf(degName, "%s.deg", txtName);
	
	inName = (char*)malloc(1+strlen(txtName) + 50);
	sprintf(inName, "%s.labelin", txtName);
	outName = (char*)malloc(1+strlen(txtName) + 50);
	sprintf(outName, "%s.labelout", txtName);

	timer tm;
	
	loadIndex();
	
	printf("load time %lf (ms)\n", tm.getTime()*1000); fflush(stdout);
	
	//memset(dagScore, 0, sizeof(dagScore));
	//memset(minPathScore, 0, sizeof(minPathScore));
	memset(dagTimeTaken,0,sizeof(dagTimeTaken));
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

	marked.resize(nodeID.size());
	int count_node = 0;
	for( int i =0 ; i < nodeID.size() ; i++) {
		if( nodeToLatLon[ nodeID[ i ] ].first < 40.6097  ) {
			marked[ i ] = 1;
			count_node++;
		}
		else 
			marked[ i ] = 0;
	}
	printf("Number of marked node: %d\n",count_node);

	frequentPickup.resize( nodeID.size(), false);
	frequentDrop.resize( nodeID.size(), false);
	
	// intializeFrequent();
	// cout<<"SIZE OF FREQUENT NODES"<<distanceFrequentNodes.size()<<endl;


	vector< vector< pair<double, vector<long long int> > > > extendEdge(nodeID.size()); 
	for( int i = 0; i < nodeID.size(); i++) {
		vector<long long int> path;
		vector< pair< double, vector<long long int> > > paths;
		extendEdges(i, i, path, paths, 0);
		extendEdge[ i ] = paths;
	}
	printf("Extend edge assigned\n"); fflush(stdout);

	//srand (time(NULL));
	srand(0);
	// 40.6097, 40.6475516, 40.6871589

	tripNumber = 0;
	for( int i=0; tripNumber < TOTAL_TRIP; i += 1) {
		cout<<tripNumber<<endl;
		long long int source, timeS, destination;
		
		do {
			source = rand() % nodeID.size() ;
		} while( !marked[ source ] );  
		
		if( !sourceTimeDestination[ source ].size() )
			continue; 

		timeS = rand()%sourceTimeDestination[ source ].size() ;

		map<long long int, vector<long long int> >::iterator it = sourceTimeDestination[ source ].begin();
		while( timeS-- ) 
			it++;

		timeS = it->first;
		
		bool found = 0;
		for( int j = 0; j < sourceTimeDestination[ source ][ timeS ].size(); j++) {
			destination = sourceTimeDestination[ source ][ timeS ][ j ];
			if( marked[ destination ] ) {
				found = 1;
				break; 
			}
		}

		if( !found )
			continue;

		bool reached = potentialScore(source, destination, timeS, extendEdge);
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

