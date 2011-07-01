
#ifndef __pvec_h
#define __pvec_h
// #define BOOST_UBLAS_NDEBUG
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <EigenvalueDecomposition.hpp>
#include <cassert>
#include <algorithm>

//#define USE_CBLAS
#ifdef USE_CBLAS
extern "C" {
#include <cblas.h>
}
#endif


#include "parsimony.h"
#include "ivymike/stupid_ptr.h"
#include "ivymike/algorithm.h"
class pvec_cgap {
    //     aligned_buffer<parsimony_state> v;
    std::vector<parsimony_state> v;
    std::vector<int> auxv;

    
    
public:
    void init( const std::vector<uint8_t> &seq ) {
//         assert( v.size() == 0 );
        v.resize(seq.size());
        auxv.resize( seq.size() );
        std::transform( seq.begin(), seq.end(), v.begin(), dna_parsimony_mapping::d2p );
        std::transform( seq.begin(), seq.end(), auxv.begin(), dna_parsimony_mapping::d2aux );
    }

    static void newview( pvec_cgap &p, pvec_cgap &c1, pvec_cgap &c2, double /*z1*/, double /*z2*/, tip_case tc ) {
        assert( c1.v.size() == c2.v.size() );

        if( c1.v.size() != c2.v.size() ) {
            throw std::runtime_error( "newview: vectors have different lengths (illegal incremetal newview on modified data?)" );
        }
        
//         p.v.resize(0);
        p.v.resize(c1.v.size());
        p.auxv.resize(c1.auxv.size());

        

        for( size_t i = 0; i < c1.v.size(); i++ ) {
            parsimony_state ps = c1.v[i] & c2.v[i];

            if( ps == 0 ) {
                ps = c1.v[i] | c2.v[i];
            }

            //p.v.push_back( ps );
            p.v[i] = ps;



            const int a1 = c1.auxv[i];
            const int a2 = c2.auxv[i];

            const bool cgap1 = (a1 & AUX_CGAP) != 0;
            const bool cgap2 = (a2 & AUX_CGAP) != 0;

//             const bool open1 = (a1 & AUX_OPEN) != 0;
//             const bool open2 = (a2 & AUX_OPEN) != 0;

            p.auxv[i] = 0;

            if( tc == TIP_TIP ) {
                if( cgap1 && cgap2 ) {
                    p.auxv[i] = AUX_CGAP;
                } else if( cgap1 != cgap2 ) {
                    p.auxv[i] = AUX_CGAP | AUX_OPEN;
                }
            } else if( tc == TIP_INNER ) {
                if( cgap1 && cgap2 ) {
                    p.auxv[i] = AUX_CGAP;
                } else if( cgap1 != cgap2 ) {
                    p.auxv[i] = AUX_CGAP | AUX_OPEN;
                }
            } else {
                if( a1 == AUX_CGAP && a2 == AUX_CGAP ) {
                    p.auxv[i] = AUX_CGAP;
                } else if( a1 == AUX_CGAP || a2 == AUX_CGAP ) {
                    p.auxv[i] = AUX_CGAP | AUX_OPEN;
                }
            }

        }
    }

    inline size_t size() {
        return v.size();
    }

    template<typename vt>
    inline void to_int_vec( std::vector<vt> &outv ) {

        outv.resize( v.size() );

        std::copy( v.begin(), v.end(), outv.begin() );
    }
    template<typename oiter_, size_t STRIDE>
    inline void to_int_vec_strided( oiter_ out ) {
        //outv.assign( v.begin(), v.end() );
        for( std::vector< parsimony_state >::iterator it = v.begin(); it != v.end(); ++it, out += STRIDE ) {
            *out = *it;
            
        }
        
//         outv.resize( v.size() );
//         std::copy( v.begin(), v.end(), outv.begin() );

    }
    
    
    template<typename vt>
    inline void to_aux_vec( std::vector<vt> &outv ) {
        //         std::cout << "v: " << v.size() << "\n";
        
        outv.resize( v.size() );
        std::copy( auxv.begin(), auxv.end(), outv.begin() );
        

//         std::for_each( auxv.begin(), auxv.end(), ostream_test(std::cout) );


    }
    
    template<typename oiter_, size_t STRIDE>
    inline void to_aux_vec_strided( oiter_ out ) {
        //         std::cout << "v: " << v.size() << "\n";
        
        for( std::vector< int >::iterator it = auxv.begin(); it != auxv.end(); ++it, out += STRIDE ) {
            if( *it == AUX_CGAP ) {
                *out = 0xFFFF;
            } else {
                *out = 0x0;
            }
        }
       // std::transform( auxv.begin(), auxv.end(), out );
        
        
        

//         std::for_each( auxv.begin(), auxv.end(), ostream_test(std::cout) );


    }
};

using namespace boost::numeric;

class probgap_model {
    ublas::matrix<double> m_evecs;
    ublas::vector<double> m_evals;
    ublas::matrix<double> m_evecs_inv;

    ublas::diagonal_matrix<double> m_evals_diag; // used temporarily.
    ublas::matrix<double> m_prob_matrix;
    ublas::matrix<double> m_temp_matrix;
    double calc_gap_freq ( const std::vector< std::vector< uint8_t > > &seqs ) {
        size_t ngaps = 0;
        size_t nres = 0;

        for( std::vector< std::vector< uint8_t > >::const_iterator it = seqs.begin(); it != seqs.end(); ++it ) {
            nres += it->size();
            ngaps += std::count_if( it->begin(), it->end(), dna_parsimony_mapping::is_gap );
        }

        double rgap = double(ngaps) / nres;
        std::cout << "gap rate: " << ngaps << " " << nres << "\n";
        std::cout << "gap rate: " << rgap << "\n";
        return rgap;
    }

public:
    probgap_model( const std::vector< std::vector<uint8_t> > &seqs ) : m_evals_diag(2), m_prob_matrix(2,2), m_temp_matrix(2,2) {
        // initialize probgap model from input sequences

        double gap_freq = calc_gap_freq( seqs );

        double f[2] = {1-gap_freq, gap_freq};

        ublas::matrix<double> rate_matrix(2,2);
        rate_matrix(0,0) = -f[0];
        rate_matrix(0,1) = f[0];
        rate_matrix(1,0) = f[1];
        rate_matrix(1,1) = -f[1];

        ublas::EigenvalueDecomposition ed(rate_matrix);

        m_evecs = ed.getV();
        m_evals = ed.getRealEigenvalues();

        // use builtin ublas lu-factorization rather than jama
        {
            ublas::matrix<double> A(m_evecs);
            ublas::permutation_matrix<size_t> pm( A.size1() );
            size_t res = ublas::lu_factorize(A,pm);
            if( res != 0 ) {
                throw std::runtime_error( " ublas::lu_factorize failed" );
            }
            m_evecs_inv = ublas::identity_matrix<double>( A.size1());

            ublas::lu_substitute(A,pm,m_evecs_inv);
        }
//         ublas::LUDecomposition lud(m_evecs);
//         m_evecs_inv = lud.pseudoinverse();
//         std::cout << "inv jama: " << m_evecs_inv << "\n";
    }

    const ublas::matrix<double> &setup_pmatrix( double t ) {

        for( int i = 0; i < 2; i++ ) {
            m_evals_diag(i,i) = exp( t * m_evals(i));
        }

//         m_prob_matrix = ublas::prod( m_evecs, m_evals_diag );
//         m_prob_matrix = ublas::prod( m_prob_matrix, m_evecs_inv );
        m_prob_matrix = ublas::prod( ublas::prod( m_evecs, m_evals_diag, m_temp_matrix ), m_evecs_inv );


//         std::cout << "pmatrix: " << m_prob_matrix << "\n";



//         throw std::runtime_error( "xxx" );

        return m_prob_matrix;
    }




};

class pvec_pgap {
    //     aligned_buffer<parsimony_state> v;
    std::vector<parsimony_state> v;
    ublas::matrix<double> gap_prob;

    
    
public:
    static ivy_mike::stupid_ptr<probgap_model> pgap_model;

    void init( const std::vector<uint8_t> &seq ) {
        assert( v.size() == 0 );
        v.resize(seq.size());

        std::transform( seq.begin(), seq.end(), v.begin(), dna_parsimony_mapping::d2p );
        //std::transform( seq.begin(), seq.end(), auxv.begin(), dna_parsimony_mapping::d2aux );

        gap_prob.resize(2, seq.size());

        for( size_t i = 0; i < seq.size(); i++ ) {
            // ( 0 )
            // ( 1 ) means gap


            if( v[i] == 0xf ) {
                gap_prob( 0, i ) = 0.0;
                gap_prob( 1, i ) = 1.0;
            } else {
                gap_prob( 0, i ) = 1.0;
                gap_prob( 1, i ) = 0.0;
            }
        }

    }

    static void newview( pvec_pgap &p, pvec_pgap &c1, pvec_pgap &c2, double z1, double z2, tip_case tc ) {
        assert( c1.v.size() == c2.v.size() );

//         p.v.resize(0);
        p.v.resize(c1.v.size());
        //p.gap_prob.resize(2, c1.v.size());

        assert( pgap_model.is_valid_ptr() );

        ublas::matrix<double> p1 = pgap_model->setup_pmatrix(z1);
        ublas::matrix<double> p2 = pgap_model->setup_pmatrix(z2);
#if 1
        ublas::matrix<double> t1 = ublas::prod(p1, c1.gap_prob);
        ublas::matrix<double> t2 = ublas::prod(p2, c2.gap_prob);
#else
        ublas::matrix<double> t1(c1.gap_prob.size1(), c1.gap_prob.size2(),0 );
        ublas::matrix<double> t2(c2.gap_prob.size1(), c2.gap_prob.size2(),0 );    
        
        {
            double *p1x = p1.data().begin();
            double *p2x = p2.data().begin();
            
            double *c1x = c1.gap_prob.data().begin();
            double *c2x = c2.gap_prob.data().begin();
            
            double *t1bx = t1.data().begin();
            double *t2bx = t2.data().begin();
            
            size_t m = p1.size2();
            size_t n = c1.gap_prob.size2();
            size_t k = p1.size1();
            
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1, p1x, m, c1x, n, 1, t1bx, n);
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1, p2x, m, c2x, n, 1, t2bx, n);
        }
#endif   
        p.gap_prob = ublas::element_prod( t1, t2 );
        //ublas::matrix< double >::const_iterator1 xxx = p.gap_prob.begin1();
       // xxx.

        //std::cout << "pvec: " << *xxx << "\n";

//         ublas::matrix< double >::iterator1 tit1 = p.gap_prob.begin1();
//
//
//         boost::io::ios_all_saver ioss(std::cout);
//         std::cout << std::setprecision(2);
//
//         std::cout << std::left;
//
// //         std::transform( tit1.begin(), tit1.end(), tit1.begin(), log );
//         std::copy( tit1.begin(), tit1.end(), std::ostream_iterator<double>(std::cout, " "));
//         std::cout << "\n";
//         ++tit1;
// //         std::transform( tit1.begin(), tit1.end(), tit1.begin(), log );
//         std::copy( tit1.begin(), tit1.end() + 4, std::ostream_iterator<double>(std::cout, " "));
//         std::cout << "\n";


        for( size_t i = 0; i < c1.v.size(); i++ ) {
            parsimony_state ps = c1.v[i] & c2.v[i];

            if( ps == 0 ) {
                ps = c1.v[i] | c2.v[i];
            }

            //p.v.push_back( ps );
            p.v[i] = ps;
        }
    }

    inline size_t size() {
        return v.size();
    }
    
    
    
    template<typename vt>
    inline void to_int_vec( std::vector<vt> &outv ) {
        outv.assign( v.begin(), v.end() );
//         outv.resize( v.size() );
//         std::copy( v.begin(), v.end(), outv.begin() );

    }
    
    
   

    template<typename vt>
    inline void to_aux_vec( std::vector<vt> &outv ) {
    
//         std::cout << "v: " << v.size() << "\n";

        outv.clear();
        outv.reserve(v.size());
//         outv.resize( v.size() );
//         std::fill( outv.begin(), outv.end(), 0 );

//         std::for_each( auxv.begin(), auxv.end(), ostream_test(std::cout) );



        ublas::matrix<double> t = get_pgap();

        // yeah! metaprogramming massacre!

        ublas::matrix< double >::iterator1 tit1 = t.begin1();
#if 0
        std::vector<double> odds;
        odds.reserve(t.size2());

//         NOTE: t.begin1() + 1).begin() is the boost::ublas way of saying "iterator over the second row"
        ivy_mike::binary_twizzle( t.begin2(), t.end2(), (t.begin1() + 1).begin(), std::back_inserter(odds), std::divides<double>() );
        
        std::transform( odds.begin(), odds.end(), odds.begin(), std::ptr_fun<double>(log) );

        //std::transform( odds.begin(), odds.end(), std::back_inserter(outv), std::bind1st( std::greater<double>(), -1.0 ) );
        
        const double odds_threshold = 0.0;
        std::transform( odds.begin(), odds.end(), std::back_inserter(outv), std::bind2nd( std::less<double>(), odds_threshold ) );
        //         std::fill( outv.begin(), outv.end(), 1 );
#else
        // set outv[i] == 1 if, prob(non-gap) is less than prob(gap)
        // this is equivalent to the above code with odds_threshold == 0.0
        ivy_mike::binary_twizzle( t.begin2(), t.end2(), (t.begin1() + 1).begin(), std::back_inserter(outv), std::less<double>() );
        
                
#endif


#if 0
//         std::transform( t.begin2(), t.end2(), odds.begin(), calc_odds<ublas::matrix< double > > );

        //boost::io::ios_all_saver ioss(std::cout);
        std::cout << std::setprecision(0);
//         std::cout << std::setw(5);
        std::cout << std::fixed;
        std::cout << std::left;

//         std::cout << t(0,0) << t(0,1) << t(0,2) <<t(0,3) <<t(0,4) << "\n";


        //std::copy( odds.begin(), odds.end(), std::ostream_iterator<double>(std::cout, " "));
        std::copy( bias.begin(), bias.end(), std::ostream_iterator<int>(std::cout));
        std::cout << "\n";

        std::copy( lomag.begin(), lomag.end(), std::ostream_iterator<unsigned char>(std::cout));
        std::cout << "\n";



        std::transform( tit1.begin(), tit1.end(), tit1.begin(), log );
        std::copy( tit1.begin(), tit1.end(), std::ostream_iterator<double>(std::cout, " "));
        std::cout << "\n";
        ++tit1;
        std::transform( tit1.begin(), tit1.end(), tit1.begin(), log );
        std::copy( tit1.begin(), tit1.end(), std::ostream_iterator<double>(std::cout, " "));
        std::cout << "\n";
#endif

    }
//     template<typename oiter_, size_t STRIDE>
//     class strided_inserter {
//         typedef std::iterator_traits<oiter_> traits;
//         typedef typename traits::value_type reference_type;
//         
//         oiter_ out;
//     public:
//         strided_inserter( oiter_ out_ ) : out(out_) {}
//         
//         strided_inserter<oiter_,STRIDE> &operator=( reference_type v ) {
//             *out = v;
//         }
//         strided_inserter<oiter_,STRIDE> &operator*() {
//             return *this;
//         }
//         
//         strided_inserter<oiter_,STRIDE> &operator++() {
//             out+=STRIDE;
//             return *this;
//         }
//         
//         strided_inserter<oiter_,STRIDE> operator++(int w) {
//             out+=STRIDE * w;
//             return *this;
//         }
//         
//         
//     };
//     
//     template<typename oiter_, size_t STRIDE>
//     inline void to_aux_vec( oiter_ out ) {
// 
//         ublas::matrix<double> t = get_pgap();
// 
//         // yeah! metaprogramming massacre!
// 
//         ublas::matrix< double >::iterator1 tit1 = t.begin1();
// 
//         // set outv[i] == 1 if, prob(non-gap) is less than prob(gap)
//         // this is equivalent to the above code with odds_threshold == 0.0
//         ivy_mike::binary_twizzle( t.begin2(), t.end2(), (t.begin1() + 1).begin(), strided_inserter<oiter_,STRIDE>(out), std::less<double>() );
//     }
    
    const ublas::matrix<double> get_pgap() {
        return gap_prob;
    }


};
ivy_mike::stupid_ptr<probgap_model> pvec_pgap::pgap_model;
#endif