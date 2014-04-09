#ifndef SIMPLEPWALIGN_HPP
#define SIMPLEPWALIGN_HPP

#include <string>

namespace simplepwalign
{

struct PairwiseAlignment
{
    std::string aligned_ref;
    std::string aligned_qry;
    int score;
};

PairwiseAlignment pairwise_global(const std::string& ref, const std::string& qry,
                                  const int match_score = 1, const int mismatch_score = -4,
                                  const int gap_open = -6, const int gap_extend = -1,
                                  const bool boundary_gaps_qry=false);


}

#endif
