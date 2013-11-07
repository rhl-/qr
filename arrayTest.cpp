#include <iostream>
#include <vector>
#include "matrixops.h"
using namespace std;
int main ()
{
vector< vector<int> > a; 
vector<int> p; 
p.push_back(1); 
p.push_back(2); 
p.push_back(3); 
vector<int> q; 
q.push_back(4); 
q.push_back(5); 
q.push_back(6); 
 
a.push_back(p); 
a.push_back(q); 

vector<int> b;
b.push_back(3);
b.push_back(2);
b.push_back(1);

vector<int> c;
matVecMult(a, b, c);
for (int i = 0; i < a.size(); i++)
{
for (int j = 0; j < a[0].size(); j++)
{
cout << a[i][j] << " ";
}
cout << endl;
}

for (vector< vector<int> >::size_type u = 0; u < c.size(); u++) 
{  
cout << c[u] << " "; 
} 
cout << endl; 

cout << c.size() << endl;
 
return 0;
}

