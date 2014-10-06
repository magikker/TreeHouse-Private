#ifndef GULO_GRAPH_WEAKLYCONNECTEDCOMPONENTSUNDIRECTED_HPP
#define GULO_GRAPH_WEAKLYCONNECTEDCOMPONENTSUNDIRECTED_HPP

#include "DepthFirst"
#include <vector>
#include <set>

namespace Gulo
{
    namespace Detail
    {
        template <typename VertexPointer>
        struct UndirectedWeaklyConnectedComponentGatherer : public DepthFirstVisitor
        {
            void startTraversal( VertexPointer v) 
            {
                components.push_back({});
            }

            void startVertex(VertexPointer v) 
            {
                components.back().insert(v);
            }

            std::vector<std::set<VertexPointer>> components;
        };

        template <typename VertexPointer>
        struct UndirectedWeaklyConnectedGraphTester : public DepthFirstVisitor
        {
            void startTraversal( VertexPointer v) 
            {
                ++started;
            }

            bool terminate() {
                return started > 1;
            }

            int started{0};
        };
    }
}

#endif // GULO_GRAPH_WEAKLYCONNECTEDCOMPONENTSUNDIRECTED_HPP
