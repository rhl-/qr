#ifndef T10_MATRIX_OPERATIONS
#define T10_MATRIX_OPERATIONS
#include <iostream>
namespace t10 {
//should be above
//template< ... >
//void matVecMult. also camelcase is for chumps and is ugly. use mat_vec_mult(A,b,c).
//also variable names matter, who the fuck writes a*b = c. Its c = A*b; so stick to  standards.
//GOOD: use of const, templates, etc.
//it doesn't really matter but you might want to use typename instead of class.
template<class Matrix, class Vector> void matVecMult (const Matrix & a, const Vector & b, Vector & c) 
{
	//who said size returns an int. -3 is an int. no vector has size -3
	//the size of a container in C is given by size_t, in C++ std::size_t
	//usually std::size_t is just unsigned int.
	int m = a.size(); //further these values are never modified so you should 
	int n = b.size(); //declare them as such, by const std::size_t a.size();
	//more generically you would write:
	// typedef typename Vector::size_t Vector_size; 
	// Vector_size m = a.size(); 
	// and you would assume vectors ``export" the size_t type.
	// type exporting is a theme.
	
	int temp; // So much for correct code, i didn't know we were ignoring the decimal portions of our results
	//However you shouldn't assume the result is a double. Assume instead that Vector exports a value_type, and use that.
	
        //also it is bad form to stick this int here. stick it inside the next for loop when you initialize it to zero.
	//this is *not* a performance hit, since no malloc is called here. this data is from the stack and declaration placement
	//will not affect runtime other than by affecting correctness. 
	for (int i = 0; i < m; i++) 
	{
		temp = 0; // all the above paragraph just to say write "Number_type temp = 0;" on this line
		for (int j = 0; j < n; j++)
		{
			temp += a[i][j]*b[j];
		}
		c.push_back(temp) //why are you pushing back here? this is really inefficient. we will talk about this in person.;
	}
}	
	/*
 *	//Rewrite version 1:
 * 	template<class Matrix, class Vector> 
	void mat_vec_mult (const Matrix & A, const Vector & b, Vector & c) {
 * 		typedef typename Vector::value_type Number;
 *		for (std::size_t i = 0; i < c.size(); ++i){
 *			c[i] = 0;
 *			for( std::size_t j = 0; j < b.size(); ++j){
 *				temp += A[i][j]*b[j];
 *			}
 *		}
 *	}
 *
 * */
} //end namespace T10 namespace
#endif //ifndef T10_MATRIX_OPERATIONS
