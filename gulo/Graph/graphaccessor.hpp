#ifndef GULO_GRAPH_GRAPHACCESSOR_HPP
#define GULO_GRAPH_GRAPHACCESSOR_HPP

namespace Gulo
{
    namespace Detail
    {
        class GraphAccessor
        {
        protected:

            template <typename Graph>
            void increaseGraphEdgeCount(Graph * graph)
            {
                ++(graph->m_edgeCount);
            }

            template <typename Graph>
            void decreaseGraphEdgeCount(Graph * graph)
            {
                --(graph->m_edgeCount);
            }

            template <typename Graph>
            void increaseGraphVertexCount(Graph * graph)
            {
                ++(graph->m_vertexCount);
            }

            template <typename Graph>
            void decreaseGraphVertexCount(Graph * graph)
            {
                --(graph->m_vertexCount);
            }

            template <typename Graph>
            static typename Graph::MemoryPool * edgePool(Graph * graph)
            {
                return graph->m_edgePool;
            }

            template <typename Graph>
            static typename Graph::MemoryPool * vertexPool(Graph * graph)
            {
                return graph->m_vertexPool;
            }
        };
    }
}

#endif // GULO_GRAPH_GRAPHACCESSOR_HPP