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
            return a[0] < b[0];
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
			temp_pair[i][0]=(4-i)/2;
			temp_pair[i][1]=2*i;
			temp_pair[i][2]=3*i;
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
		int size=4;
		int i=0;
	 	do
		{
			if (temp_pair[i][0]==temp_pair[i+1][0])
				{
					temp_pair.erase( temp_pair.begin() + i);
					i-=1;
					size-=1;
				}
			i=i+1;			
		} while (i<size-1);
		
			 	cout<<"After Unique "<<endl;

	 		for (int i=0; i<size; i++) 
		{
					cout<<temp_pair[i][0]<<'\t'<<temp_pair[i][1]<<'\t'<<temp_pair[i][2]<<'\t'<<endl;
		}
	return 0;
}
