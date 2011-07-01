#include <stdint.h>
#include <iostream>
// #include "ivymike/sdf.h"

class test_impl1 {
public:
    struct xxx {
        int a;
        int b;
        xxx( int a_ ) : a(a_), b(0) {}
        
        void print() {
            std::cout << "impl1: " << a << " " << b << "\n";
        }
    };
    
};



class test_impl2 {
public:
    struct xxx {
        int8_t a;
        
        xxx( int a_ ) : a(a_) {}
        
        void print() {
            std::cout << "impl2: " << int(a) << "\n";
        }
    };
    
};



template<typename impl_>
class test : private impl_ {
public:
    using typename impl_::xxx;
    
};

typedef test<test_impl1> test1;
typedef test<test_impl2> test2;


int main() {
    test1 x;
    test1::xxx y(129);
    
    test2::xxx z(129);
    y.print();
    z.print();
    
//     ivy_mike::sdf_full sdf( std::cin, true );
    
}