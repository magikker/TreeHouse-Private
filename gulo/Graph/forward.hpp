#ifndef GULO_GRAPH_FORWARD_HPP
#define GULO_GRAPH_FORWARD_HPP

namespace Gulo
{
    namespace Detail
    {
        template <typename Edge, typename DirectedOption>
        class VertexExtension {}; 

        template <typename Vertex, typename Edge, typename DirectedOption>
        class EdgeExtension {};

        template <typename Graph, typename Vertex, typename DirectedOption>
        class GraphExtension {};

        template <typename DirectedOption>
        struct ModifyVertexStorage;

        template <typename DirectedOption>
        struct CopyEdgeExtension;

        template <typename DirectedOption>
        struct CopyVertexExtension;

        template <typename DirectedOption>
        struct GraphContraction;

        template <typename Directionality>
        struct TraverseOutgoingEdges;

        template <typename DirectedOption>
        class GraphTraversal;

        template <typename Edge, typename DirectedOption>
        struct ParallelEdgeDetection;

        template <typename DirectedOption>
        struct StandardEdgeIsomorphismTest;

        template <typename DerivedGraph, typename DerivedVertex, typename DerivedEdge, typename DirectedOption>
        class Graph;

        template <typename DataT, typename CustomDataT, typename DerivedGraph, typename DerivedVertex, typename DerivedEdge, typename DirectedOption>
        class GraphVertex;

        template <typename DataT, typename CustomDataT, typename DerivedGraph, typename DerivedVertex, typename DerivedEdge, typename DirectedOption>
        class GraphEdge;

        template <typename T> 
        struct InOutUndirected;
    }

    struct StoreOutEdges;
    
    struct StoreInOutEdges;
}        

#endif // GULO_GRAPH_FORWARD_HPP
