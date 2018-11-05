# SHARE
This project can be run on any city map and we have experimented and testes on the data for the cities: Singapore, New York and San Francisco

Steps to run the project:
1. Download the Data Set for the city for which you want to run the simulation from the following link: 
    https://drive.google.com/drive/folders/1DiVSOqANI3Ww0jHSKMw5VcxH04aIJsHC?usp=sharing

2. Keep all the data files in the same folder as that of all the other code.

3. To compile the Project run: 
    " g++ -O3 -std=c++11 realtimesimulation.cpp Hungarian.cpp "
  
    This will generate the 'a.out' file and then execute it with parameters.
    
4. To run :
    Usage: ./a.out [location={SF|NY|SG}] [alpha=1.3] [route={dij|dag|dex}] [maxDepth=0.2] [assign={hun|pxa}] [maxCab=2000] [cabCapacity={2|3}] [Optional: starttime endtime]
    
    The parameters can be explained as follows:
    Location: SF|NY|SG represents the city in which you are running the simulation. SF= San Francisco, NY= New York , SG= Singapore 
    
    Alpha= Percentage Amount of extra distance that we allow the passenger to travel as compared to the shortest path from its            source to its destination.
    
    Route: This parameter chooses the type of route recommendation algorithm that we run. 
    
            'dij' = Dijkstra (Shortest Path)
            
            'dag' = Use the Dynamic Programming algorithm after converting the graph into a DAG.
            
            'dex' = Use the Dynamic Programming algorithm after converting the graph into a DAG but along with back edges of                length equal to Max Depth.
    
    Max Depth: It is the maximum length of the reverse(back) edges allowed in the DAG.
    
    Assign: This represents the assignment algorithm that is used for assignment of free cabs to passengers.
    
            'hun' = Hungarian Algorithm for empty cab assignment
            
            'pxa' = Greedy Price assignment algorithm used for comparison with the VLDB paper.
    
    MaxCab = Maximum number of cabs running in the city
    
    cabCapacity = This represents the maximum number of passengers in the cab.
    
    Optional: starttime -> The start time of the simulation in the day
    
              endtime -> The end time of the simulation in the day
    
    
    Thus, a sample command would look like: ./a.out NY 1.2 dex 0.3 hun 4000 3 480 570
    Here, the definitions for each variable have been mentioned in the paper.
    
    
