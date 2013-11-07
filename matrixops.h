#include <iostream>

template<class Matrix, class Vector> void matVecMult (const Matrix & a, const Vector & b, Vector & c) 
{
	int m = a.size();
	int n = b.size();
	int temp;
	for (int i = 0; i < m; i++)
	{
		temp = 0;
		for (int j = 0; j < n; j++)
		{
			temp += a[i][j]*b[j];
		}
		c.push_back(temp);
	}
}	
