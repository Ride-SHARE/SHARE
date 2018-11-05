#include <bits/stdc++.h>

using namespace std;


char * txtName;
string location;
int locationFactor = 1;

int n;

int INF = (1<<30);

int main(int argc, char const *argv	[])
{
	if (argc < 3) {
		printf("Usage: ./a.out [location={BJ|SF|NY}] [connected=0/strongly connected=1]\n");
		return 0;
	}
	location = argv[1];
	int opt = 0;
	sscanf(argv[2],"%d",&opt);
	txtName = (char*)malloc(50);
	if (location.compare("BJ") == 0) {
		sprintf(txtName, "beijingIndex");
	}
	else if (location.compare("NY") == 0) {
		sprintf(txtName, "nyIndex");
		locationFactor = 2;
	}
	else if (location.compare("SF") == 0) {
		sprintf(txtName, "sfIndex");	}
	else {
		printf("Location not recognized.");
		return 1;
	}

	ifstream file;
	file.open( txtName );
	string s;
	getline(file, s);
	stringstream ss( s );
	ss >> n;
	vector<vector<pair<int,int> > > e(n);

	while( getline(file, s) ) {
		stringstream ss( s );
		int a,b,w;
		ss>> a >> b >> w;
		//printf("%d %d %d\n",a,b,w);
		e[a].push_back(make_pair(b,w));
		if (opt == 0) {
			e[b].push_back(make_pair(a,w));
		}
	}

	for (int i=0; i<n; i++) {
		vector<int> dp(n);
		priority_queue<pair<int,int> > pq;

		for (int j=0; j<n; j++) {
			dp[j] = INF;
		}
		dp[i] = 0;
		pq.push(make_pair(0,i));

		while (!pq.empty()) {
			int at = pq.top().second;
			int w = -pq.top().first;
			//printf(" At %d %d\n",at,w);
			pq.pop();
			for (int i=0; i<e[at].size(); i++) {
				if (w + e[at][i].second < dp[e[at][i].first]) {
					dp[e[at][i].first] = w + e[at][i].second;
					//printf("  Added %d %d\n",e[at][i].first,dp[e[at][i].first]);
					pq.push(make_pair(-dp[e[at][i].first], e[at][i].first));
				}
			}
		}

		printf("Not reachable from %d:",i);
		for (int j=0; j<n; j++) {
			if (dp[j] == INF) {
				printf(" %d",j);
			}
		}
		printf("\n");
	}
	
	return 0;
}

