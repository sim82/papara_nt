
#include <iostream>
#include <iterator>
#include <tbb/task.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/banded.hpp>

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

namespace {
using namespace boost::numeric::ublas;

static void mytred2( matrix<double> &a, const int n, vector<double> &d, vector<double> &e)
{
    int     l, k, j, i;
    double  scale, hh, h, g, f;

    for (i = n; i > 1; i--)
    {
        l = i - 1;
        h = 0.0;
        scale = 0.0;

        if (l > 1)
        {
            for (k = 1; k <= l; k++)
                scale += fabs(a(k - 1, i - 1));
            if (scale == 0.0)
                e[i - 1] = a(l - 1, i - 1);
            else
            {
                for (k = 1; k <= l; k++)
                {
                    a(k - 1, i - 1) /= scale;
                    h += a(k - 1, i - 1) * a(k - 1, i - 1);
                }
                f = a(l - 1,i - 1);
                g = ((f > 0) ? -sqrt(h) : sqrt(h)); /* diff */
                e[i - 1] = scale * g;
                h -= f * g;
                a(l - 1, i - 1) = f - g;
                f = 0.0;
                for (j = 1; j <= l; j++)
                {
                    a(i - 1,j - 1) = a(j - 1,i - 1) / h;
                    g = 0.0;
                    for (k = 1; k <= j; k++)
                        g += a(k - 1, j - 1) * a(k - 1, i - 1);
                    for (k = j + 1; k <= l; k++)
                        g += a(j - 1,k - 1) * a(k - 1,i - 1);
                    e[j - 1] = g / h;
                    f += e[j - 1] * a(j - 1, i - 1);
                }
                hh = f / (h + h);
                for (j = 1; j <= l; j++)
                {
                    f = a(j - 1, i - 1);
                    g = e[j - 1] - hh * f;
                    e[j - 1] = g;
                    for (k = 1; k <= j; k++)
                        a(k - 1, j - 1) -= (f * e[k - 1] + g * a(k - 1, i - 1));
                }
            }
        }
        else
            e[i - 1] = a(l - 1, i - 1);
        d[i - 1] = h;
    }
    d[0] = 0.0;
    e[0] = 0.0;

    for (i = 1; i <= n; i++)
    {
        l = i - 1;
        if (d[i - 1] != 0.0)
        {
            for (j = 1; j <= l; j++)
            {
                g = 0.0;
                for (k = 1; k <= l; k++)
                    g += a(k - 1, i - 1) * a(j - 1, k - 1);
                for (k = 1; k <= l; k++)
                    a(j - 1, k - 1) -= g * a(i - 1, k - 1);
            }
        }
        d[i - 1] = a(i - 1, i - 1);
        a(i - 1, i - 1) = 1.0;
        for (j = 1; j <= l; j++)
            a(i - 1, j - 1) = a(j - 1, i - 1) = 0.0;
    }


}


// this seems to be the function tqli from NR chaper 11, but with rows/columns of z switched,
// so that the output eigenvectors are in the rows of z
static int mytqli( vector<double> &d, vector<double> &e, const int n, matrix<double> &z)
{
    int     m, l, iter, i, k;
    double  s, r, p, g, f, dd, c, b;

    for (i = 2; i <= n; i++)
        e[i - 2] = e[i - 1];

    e[n - 1] = 0.0;

    for (l = 1; l <= n; l++)
    {
        iter = 0;
        do
        {
            for (m = l; m <= n - 1; m++)
            {
                dd = fabs(d[m - 1]) + fabs(d[m]);
                if (fabs(e[m - 1]) + dd == dd)
                    break;
            }

            if (m != l)
            {
                assert(iter < 30);

                g = (d[l] - d[l - 1]) / (2.0 * e[l - 1]);
                r = sqrt((g * g) + 1.0);
                g = d[m - 1] - d[l - 1] + e[l - 1] / (g + ((g < 0)?-fabs(r):fabs(r)));/*MYSIGN(r, g));*/
                s = c = 1.0;
                p = 0.0;

                for (i = m - 1; i >= l; i--)
                {
                    f = s * e[i - 1];
                    b = c * e[i - 1];
                    if (fabs(f) >= fabs(g))
                    {
                        c = g / f;
                        r = sqrt((c * c) + 1.0);
                        e[i] = f * r;
                        c *= (s = 1.0 / r);
                    }
                    else
                    {
                        s = f / g;
                        r = sqrt((s * s) + 1.0);
                        e[i] = g * r;
                        s *= (c = 1.0 / r);
                    }
                    g = d[i] - p;
                    r = (d[i - 1] - g) * s + 2.0 * c * b;
                    p = s * r;
                    d[i] = g + p;
                    g = c * r - b;
                    for (k = 1; k <= n; k++)
                    {
                        f = z(i,k-1);
                        z(i,k-1) = s * z(i - 1, k - 1) + c * f;
                        z(i - 1, k - 1) = c * z(i - 1, k - 1) - s * f;
                    }
                }

                d[l - 1] = d[l - 1] - p;
                e[l - 1] = g;
                e[m - 1] = 0.0;
            }
        }
        while (m != l);
    }



    return (1);
}


static void makeEigen( matrix<double> &_a, const int n, vector<double> &d, vector<double> &e)
{
    mytred2(_a, n, d, e);
    mytqli(d, e, n, _a);
}
}

int main() {
    {
        long n = 10;
        long sum = 0;
        using namespace tbb;
        fib_task &a = *new(task::allocate_root()) fib_task(n, &sum);
        task::spawn_root_and_wait(a);
        std::cout << sum << "\n";
    }
    
    using namespace boost::numeric::ublas;
    
    const size_t n = 2;
    
    vector<double> frequencies(n);
    frequencies(0) = 0.55;
    frequencies(1) = 0.45;
    
    matrix<double> rate_matrix(n,n);
    
    
    for ( int i=0; i < n; i++)
    {
        for ( int j=i+1;  j < n; j++)
        {
            //             double factor =  initialRates[m++];
            double factor = 1.0; // initial rate = 1
            rate_matrix(i,j) = rate_matrix(j,i) = factor * sqrt( frequencies(i) * frequencies(j));
            rate_matrix(i,i) -= factor * frequencies(j);
            rate_matrix(j,j) -= factor * frequencies(i);
        }
    }
    
    
    vector<double> eigenvalues(2);
    vector<double> e(2);
    matrix<double> eigenvectors = rate_matrix;
    
    
    
    
    makeEigen(eigenvectors, 2, eigenvalues, e);
    
    std::cout << eigenvectors(0,0) << " " << eigenvectors(0,1) << "\n" << eigenvectors(1,0) << " " << eigenvectors(1,1) << "\n";
    diagonal_matrix<double> diag_evs(2); 
    for( int i = 0 ; i < n; i++ ) {
        diag_evs(i,i) = eigenvalues(i);
    }
    
    
    std::cout << diag_evs(0,0) << " " << diag_evs(0,1) << "\n" << diag_evs(1,0) << " " << diag_evs(1,1) << "\n";
    
//     std::copy( m[0].begin1(), m[0].end1(), std::ostream_iterator<double>( std::cout, "\n" ));
//     std::copy( m[1].begin1(), m[1].end1(), std::ostream_iterator<double>( std::cout, "\n" ));
    
//     double a = m(0,0);
//     double b = m(0,1);
//     double c = m(1,0);
//     double d = m(1,1);
//     
//     double e1 = (a+d)/2 + sqrt(4 * b * c + (a-d) * (a - d))/2;
//     double e2 = (a+d)/2 - sqrt(4 * b * c + (a-d) * (a - d))/2;
//     
//     std::copy( m.begin2(), m.end2(), std::ostream_iterator<double>( std::cout, "\n" ));
//     
//     
//     std::cout << e1 << " " << e2 << "\n";
//     
}