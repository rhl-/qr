#include <iostream>
#include <vector>
#include "matrix_operations.h"

typedef std::vector<int> Vector;
int main ()
{
std::vector< Vector > A(3,std::vector<int>(2,1));


Vector b;
b.push_back(3);
b.push_back(2);
b.push_back(1);

Vector c(2,0);
t10::mat_vec_mult(A, b, c);
for (int i = 0; i < A.size(); i++){
	for (int j = 0; j < A[0].size(); j++){
		std::cout << A[i][j] << " ";
	}
	std::cout << std::endl;
}

for (std::vector< Vector >::size_type u = 0; u < c.size(); u++) {  
	std::cout << c[u] << " " << std::endl; 
} 
std::cout << std::endl; 

std::cout << "size of c " << c.size() << std::endl;
 
return 0;
}

