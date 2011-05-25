
#include <iostream>

#include <tbb/task.h>
class fib_task : public tbb::task {
    const long m_n;
    long *const m_sum;
    
public:
    fib_task( long n, long *sum ) : m_n(n), m_sum(sum) {}
//     virtual ~fib_task() {
//         std::cout << "term\n";
//     }
    
    task * execute() {
//         std::cout << "execute\n";
        if( m_n < 2 ) {
            *m_sum = m_n;
        } else {
            long x, y;
            
            fib_task &a = *new(allocate_child() ) fib_task(m_n-1, &x);
            fib_task &b = *new(allocate_child() ) fib_task(m_n-2, &y);
            
            set_ref_count(3);
            spawn(b);
            spawn_and_wait_for_all(a);
            
            *m_sum = x + y;
            
        }
        return 0;
    }
    
};

int main() {
    long n = 42;
    long sum = 0;
    using namespace tbb;
    fib_task &a = *new(task::allocate_root()) fib_task(n, &sum);
    task::spawn_root_and_wait(a);
    std::cout << sum << "\n";
    
}