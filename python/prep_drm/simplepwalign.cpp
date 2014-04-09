#include "simplepwalign.hpp"

#include "seqan/align.h"

using namespace seqan;
using std::string;

namespace simplepwalign
{

template<typename TRow>
string aligned_string(TRow& row)
{
    typedef typename Iterator<TRow>::Type TRowIterator;

    string result(length(row), '-');

    TRowIterator it = begin(row);
    TRowIterator itEnd = end(row);
    for(int i = 0; it != itEnd; ++it, ++i) {
        if(!isGap(it))
            result[i] = *it;
    }
    return result;
}

template<bool QGAP>
PairwiseAlignment _pairwise_global(const std::string& ref, const std::string& qry,
                                   const int match_score, const int mismatch_score,
                                   const int gap_open, const int gap_extend)
{
    typedef String<char> TSequence;
    typedef Align<TSequence,ArrayGaps> TAlign;
    typedef Row<TAlign>::Type TRow;
    typedef Iterator<TRow>::Type TRowIterator;

    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align, 0), ref);
    assignSource(row(align, 1), qry);
    Score<int,Simple> score_scheme(match_score, mismatch_score, gap_extend, gap_open);
    AlignConfig<false,QGAP,false,QGAP> a_config;
    PairwiseAlignment result;
    result.score = globalAlignment(align, score_scheme, a_config);

    result.aligned_ref = aligned_string(row(align, 0));
    result.aligned_qry = aligned_string(row(align, 1));
    return result;
}


PairwiseAlignment pairwise_global(const std::string& ref, const std::string& qry,
                                  const int match_score, const int mismatch_score,
                                  const int gap_open, const int gap_extend,
                                  const bool boundary_gaps_qry)
{
    if(boundary_gaps_qry)
        return _pairwise_global<true>(ref, qry, match_score, mismatch_score,
                                      gap_open, gap_extend);
    else
        return _pairwise_global<false>(ref, qry, match_score, mismatch_score,
                                      gap_open, gap_extend);
}

}
