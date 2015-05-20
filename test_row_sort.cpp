#include <iostream>
#include <random>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <tuple>
#include <dirent.h>
#include <algorithm>
#include <functional>
#include <array>
# include "defs.h"
# include "force_structured.h"

using namespace std;

// compare two rows by elements of the first elemtn of the row for sroting 

    struct RowSort {
        bool operator()(vector<int> a, vector<int>  b)
        {   
            return a[1] < b[1];
        }   
    } ;
    
    struct RowUnique {
        bool operator()(vector<int> a, vector<int>  b)
        {   
            return a[1] != b[1];
        }   
    } ;
	
	int main() {
	vector<vector<int>> temp_pair(4,vector<int> (3))	;
	for (int i=0; i<4; i++) 
		{
			for (int j=0; j<3; j++) 
				{
					temp_pair[i][j]=(4-i)/2;
				}
		}	
		cout<<"Before Sort"<<endl;
		for (int i=0; i<4; i++) 
		{
					cout<<temp_pair[i][0]<<'\t'<<temp_pair[i][1]<<'\t'<<temp_pair[i][2]<<'\t'<<endl;
		}		
		
	 	sort (temp_pair.begin(),temp_pair.end(), RowSort());
	 	
	 	cout<<"After Sort"<<endl;

	 		for (int i=0; i<4; i++) 
		{
					cout<<temp_pair[i][0]<<'\t'<<temp_pair[i][1]<<'\t'<<temp_pair[i][2]<<'\t'<<endl;
		}
	return 0;
}
