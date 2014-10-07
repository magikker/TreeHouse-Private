#ifndef GULO_BASICGRAPH_HPP
#define GULO_BASICGRAPH_HPP

namespace Gulo
{   
    template <typename DirectedOption, typename VertexT, typename EdgeT>
    class Graph;

    template <typename DirectedOption, typename VertexT, typename EdgeT>
    class GraphEdge;

	template <typename DirectedOption, typename VertexT, typename EdgeT>
	class GraphVertex;

    template <typename DirectedOption, typename VertexT, typename EdgeT>
    class GraphVertex : public Detail::GraphVertex<VertexT, void, Graph<DirectedOption, VertexT, EdgeT>, GraphVertex<DirectedOption, VertexT, EdgeT>, GraphEdge<DirectedOption, VertexT, EdgeT>, DirectedOption >
    {
		typedef Detail::GraphVertex<VertexT, void, ::Gulo::Graph<DirectedOption, VertexT, EdgeT>, ::Gulo::GraphVertex<DirectedOption, VertexT, EdgeT>, ::Gulo::GraphEdge<DirectedOption, VertexT, EdgeT>, DirectedOption> Base;
    
    protected:
        
        using Base::Base;
        
    };

    template <typename DirectedOption, typename VertexT, typename EdgeT>
    class GraphEdge : public Detail::GraphEdge<EdgeT, void, Graph<DirectedOption, VertexT, EdgeT>, GraphVertex<DirectedOption, VertexT, EdgeT>, GraphEdge<DirectedOption, VertexT, EdgeT>, DirectedOption>
    {
		typedef Detail::GraphEdge<EdgeT, void, ::Gulo::Graph<DirectedOption, VertexT, EdgeT>, ::Gulo::GraphVertex<DirectedOption, VertexT, EdgeT>, ::Gulo::GraphEdge<DirectedOption, VertexT, EdgeT>, DirectedOption> Base;
        
    protected:

        using Base::Base;
    
    };

    template <typename DirectedOption, typename VertexT, typename EdgeT>
    class Graph: public Detail::Graph<Graph<DirectedOption, VertexT, EdgeT>, GraphVertex<DirectedOption, VertexT, EdgeT>, GraphEdge<DirectedOption, VertexT, EdgeT>, DirectedOption>
    {
		typedef Detail::Graph<::Gulo::Graph<DirectedOption, VertexT, EdgeT>, ::Gulo::GraphVertex<DirectedOption, VertexT, EdgeT>, ::Gulo::GraphEdge<DirectedOption, VertexT, EdgeT>, DirectedOption> Base;
        
    public:

        Graph()
        : Base()
        {
        }
        
        Graph(const Graph<DirectedOption, VertexT, EdgeT> & other)
        : Base(static_cast<const Base&>(other))
        {
        }
        
        Graph(Graph<DirectedOption, VertexT, EdgeT> && other)
        : Base(std::move(other))
        {
        }
        
        Graph & operator=(const Graph & other)
        {
            if (this != &other) {
                static_cast<Base&>(*this) = static_cast<const Base&>(other);
            }
            return *this;
        }
        
        Graph & operator=(Graph && other)
        {
            static_cast<Base&>(*this) = std::move(other);
            return *this;
        }
    };
}

#endif // GULO_BASICGRAPH_HPP 