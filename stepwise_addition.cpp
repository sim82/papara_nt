/*
 * Copyright (C) 2011 Simon A. Berger
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 */
#define PWDIST_INLINE
#include "pairwise_seq_distance.h"

#include "ivymike/tdmatrix.h"

#include "parsimony.h"
#include "pvec.h"
#include "fasta.h"
#include "ivymike/tree_parser.h"

using namespace ivy_mike::tree_parser_ms;
// void pairwise_seq_distance(std::vector< std::vector<uint8_t> > &seq_raw, ivy_mike::tdmatrix<int> &, const scoring_matrix &sm, const int gap_open, const int gap_extend, const int n_thread );

template<class pvec_t>
class my_adata_gen : public ivy_mike::tree_parser_ms::adata {
//     static int ct;
    //std::vector<parsimony_state> m_pvec;
    pvec_t m_pvec;

public:
//     int m_ct;
    my_adata_gen() {

//         std::cout << "my_adata\n";

    }

    virtual ~my_adata_gen() {

//         std::cout << "~my_adata\n";

    }

    virtual void visit() {
//         std::cout << "tr: " << m_ct << "\n";
    }
    void init_pvec(const std::vector< uint8_t >& seq) {


        m_pvec.init( seq );

    }
    pvec_t &get_pvec() {
        return m_pvec;
    }


};

template<class ndata_t>
class my_fact_gen : public ivy_mike::tree_parser_ms::node_data_factory {

    virtual ndata_t *alloc_adata() {

        return new ndata_t;
    }

};


class step_add {
    typedef pvec_pgap pvec_t;
    
    typedef my_adata_gen<pvec_t> my_adata;
    typedef my_fact_gen<my_adata> my_fact;
    //std::auto_ptr<ivy_mike::tree_parser_ms::ln_pool> m_ln_pool;
    ivy_mike::tree_parser_ms::ln_pool m_ln_pool;
    std::string m_seq_file_name;
    std::vector<std::string> m_qs_names;
    std::vector<std::vector<uint8_t> > m_qs_seqs;
    
    std::vector<std::vector <uint8_t> > m_qs_nongappy;
    
    ivy_mike::tdmatrix<float> m_pw_dist;
    scoring_matrix m_pw_scoring_matrix;
    
public:
    step_add( const char *seq_name ) 
    : m_ln_pool( new my_fact() ),
    m_seq_file_name(seq_name),
    m_pw_scoring_matrix(3,0)
    
    {
        {
            std::ifstream qsf( m_seq_file_name.c_str() );
            read_fasta( qsf, m_pw_scoring_matrix, m_qs_names, m_qs_seqs);
        } 
    }
    
    std::pair<size_t,size_t> calc_dist_matrix() {
        m_pw_dist.init_size(m_qs_names.size(), m_qs_names.size());
    
        
        
        ivy_mike::tdmatrix<int> out_scores(m_qs_names.size(), m_qs_names.size());
        
        pairwise_seq_distance(m_qs_seqs, out_scores, m_pw_scoring_matrix, -5, -2, 4);
     
        size_t li, lj;
        float lowest_dist = 1e8;
        
        
        for( size_t i = 0; i < out_scores.size(); i++ ) {
            
            for( size_t j = 0; j < out_scores[i].size(); j++ ) {
            
                // three modes for normalizing: min, max and mean
                //const float norm = std::min( ma[i][i], ma[j][j] );
                //             const float norm = std::max( ma[i][i], ma[j][j] );
                const float norm = (out_scores[i][i] + out_scores[j][j]) * 0.5;
                
                
                int mae;
                if( i <= j ) {
                    mae = out_scores[i][j];
                    //                 mae = ma[j][i];
                } else {
                    mae = out_scores[j][i];
                    
                }
                
                const float dist = 1.0 - (mae / norm);
                m_pw_dist[i][j] = dist;
                
                if( i != j && dist < lowest_dist ) {
                    lowest_dist = dist;
                    li = i;
                    lj = j;
                }
                
            }
        
        }
        
        return std::pair<size_t,size_t>(li,lj);
    }
    
};

int main() {
    //mapped_file qsf( "test_1604/1604.fa" );
	
    
    
    
    

    
    step_add sa("test_1604/1604.fa.100");
    std::pair<size_t,size_t> start_pair = sa.calc_dist_matrix();
    
    std::cout << "start: " << start_pair.first << " " << start_pair.second << "\n";
    
}