#ifndef __raxml_interface_h
#define __raxml_interface_h


#include "ivymike/tree_parser.h"
#include <string>
#include <vector>
#include <map>
#include <stdint.h>

void optimize_branch_lengths( ivy_mike::tree_parser_ms::lnode *tree, const std::map<std::string, const std::vector<uint8_t> * const> &name_to_seq );
#endif
