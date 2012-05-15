// just a test if i can write this faster in c++ than in ruby (I always forget half of the ruby syntax and hate to look up the details again and again...)

#include <iostream>
#include <cstdlib>

int main( int argc, char *argv[] ) {
	

	std::string line;

	
    bool fixed_number = false;
    size_t num = 0;
    
    
	double p;
	if( argc > 1 ) {
		p = atof(argv[1]);
		if( p < 0 ) {
			std::cerr << "bad p: " << p << "\n";
		} 
		
		num = size_t(p);
        fixed_number = p >= 1.0;
		
	} else {
		p = 0.1;
	}
	
	
	
	bool is_on = false;
    while( !std::cin.eof() ) {
        line.clear();
        std::getline( std::cin, line );
        
        if( line.size() > 0 ) {
            if( line[0] == '>' ) {
                is_on = (std::rand() / float(RAND_MAX)) < p;
            }
        }
        
        if( is_on ) {
            std::cout << line << "\n";
        }
        
    }
}
