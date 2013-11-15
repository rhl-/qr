#include <boost/numeric/ublas/io.hpp> //io
#include <iomanip>
#include <sstream>
namespace ublas = boost::numeric::ublas;
namespace t10 {
	template< typename Matrix>
	std::string print_matrix( const Matrix & M, std::size_t p=7){
		std::string str;
		std::stringstream ss;
		ss.setf( std::ios::fixed, std:: ios::floatfield );
		ss << std::setw( p) << std::setprecision( p) <<  M;
		str = ss.str();
		//kill size thing
		std::size_t lefts = str.find(std::string("["));
		std::size_t rights = str.find(std::string("]"));
		str.replace(lefts,rights-lefts+1,"");
		
		//start with new line	
		std::size_t found = str.find(std::string("(("));
		str.replace(found,2,std::string("\n"));
		//end with new line
		found = str.find(std::string("))"));
		str.replace(found,2,std::string("\n"));
	
		//end of each row
		found = str.find(std::string("),"));
		while( found != std::string::npos){
			str.replace(found,2,std::string("\n"));
			found = str.find(std::string("),"));
		}
		//remove beginning of row
		found = str.find(std::string("("));
		while( found != std::string::npos){
			str.replace(found,1,std::string(""));
			found = str.find(std::string("("));
		}
		//replace , with tab
		found = str.find(std::string(","));
		while( found != std::string::npos){
			str.replace(found,1,std::string("\t"));
			found = str.find(std::string(","));
		}
		
		return str;
	}



	
	template< typename Matrix>
	std::string print_matrix( const Matrix & M, bool matlab){
		std::string str;
		std::stringstream ss;
		ss << M;
		str = ss.str();
		std::size_t lefts = str.find(std::string("["));
		std::size_t rights = str.find(std::string("]"));
		str.replace(lefts,rights-lefts+1,"");
		std::size_t found = str.find(std::string("(("));
		str.replace(found,2,std::string("["));
		found = str.find(std::string("))"));
		str.replace(found,2,std::string("]"));
		found = str.find(std::string("),"));
		while( found != std::string::npos){
			str.replace(found,2,std::string(";"));
			found = str.find(std::string("),"));
		}

		found = str.find(std::string("("));
		while( found != std::string::npos){
			str.replace(found,1,std::string(""));
			found = str.find(std::string("("));
		}


		return str;	
	}
} //end namespace t10
