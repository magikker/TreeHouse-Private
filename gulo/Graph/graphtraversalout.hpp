#ifndef GULO_GRAPHTRAVERSALOUT_HPP
#define GULO_GRAPHTRAVERSALOUT_HPP

#include "graphtraversal.hpp"

namespace Gulo
{
    namespace Detail
    {
        /////////////////////////////
        // Traverse Outgoing Edges //
        /////////////////////////////

        template <template <typename> class DirectedOption>
        struct TraverseOutgoingEdges<DirectedOption<StoreOutEdges>>
        {
            template <typename Vertex, typename Func>
            static void go(Vertex* v, Func & func, const typename Vertex::Graph::Edge * prev)
            {
                for (auto * edge : v->outgoingEdges()) {
                    func(edge, edge->vertex2());
                }
            }
        };
    }
}

#endif // GULO_GRAPHTRAVERSALOUT_HPP