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
    Usage: ./a.out [location={SF|NY|SG}] [alpha=1.3] [route={dij|dag|dex}] [maxDepth=0.2] [assign={hun|gsp|pxa|sdp}] [maxCab=2000] [cabCapacity=2] [Optional: starttime endtime]
    
    Thus, a sample command would look like: ./a.out NY 1.2 dex 0.3 hun 4000 3 480 570
    Here, the definitions for each variable have been mentioned in the paper.
    
    
