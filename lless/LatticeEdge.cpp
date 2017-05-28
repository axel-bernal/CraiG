#include "LatticeEdge.h"
#include "LatticeVertex.h"

namespace lless {

    /**
     * Main Constructor
     * @param st the node to insert at this position in the lattice.
     * @param wordLength word length up to this point, from the last seen
     * synchonization state.
     * @param backPType the parsing type of the previous node
     * @param backPhase the phase of the previous node
     * @param backKth the rank of the previous node
     * @param score path score up to this point.
     */
    LatticeEdge::LatticeEdge(Tag &st,
                             int wordLength,
                             TParseNode backpType,
                             int backPhase,
                             int backKth,
                             double score) {

        this->st = st;
        this->wordLength = wordLength;
        this->_backpType = backpType;
        this->backPhase = backPhase;
        this->backKth = backKth;
        this->_score = score;
    }

    /**
     * Move back to the previous vertex in the lattice.
     */
    LatticeEdge *LatticeEdge::backtrack(LatticeVertex ***vptr) {
        return (*vptr[_backpType][backPhase])[backKth];
    }

}
