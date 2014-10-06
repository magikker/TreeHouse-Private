#ifndef GULO_GRAPH_ARBORESCENCETEST_HPP
#define GULO_GRAPH_ARBORESCENCETEST_HPP

namespace Gulo
{
    namespace Detail
    {
        template <typename VertexPointer>
        struct ArborescenceTester : public DepthFirstVisitor
        {
            void startTraversal(VertexPointer v)
            {
                if (initialVertex == nullptr) {
                    initialVertex = v;
                }
            }

            void reachVisitedVertex(VertexPointer v)
            {
                failed = true;
            }

            void reachQueuedVertex(VertexPointer v)
            {
                failed = true;
            }

            bool terminate() const
            {
                return failed;
            }

            bool failed{false};
            VertexPointer initialVertex{nullptr};
        };
    }
}

#endif // GULO_GRAPH_ARBORESCENCETEST_HPP
