





#include "fasta.h"



void pairwise_seq_distance( std::vector< std::vector<uint8_t> > &seq_raw );

int main() {
    //mapped_file qsf( "test_1604/1604.fa" );
	std::ifstream qsf( "test_1604/1604.fa" );
    std::vector<std::string> qs_names;
    std::vector<std::vector<uint8_t> > qs_seqs;
    
    std::vector<std::vector <uint8_t> > qs_nongappy;
    
    
    
    
    read_fasta( qsf, qs_names, qs_seqs);
    
    pairwise_seq_distance(qs_seqs);
    
}