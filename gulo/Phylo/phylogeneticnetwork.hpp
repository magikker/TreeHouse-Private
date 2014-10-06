#ifndef GULO_PHYLOGENETICNETWORK_HPP
#define GULO_PHYLOGENETICNETWORK_HPP

namespace Gulo
{       
    template <typename DirectedOption, typename VertexT, typename EdgeT, typename EdgeLengthT, typename VertexLabelT>
    class PhylogeneticNetworkVertex;

    template <typename DirectedOption, typename VertexT, typename EdgeT, typename EdgeLengthT, typename VertexLabelT>
    class PhylogeneticNetworkEdge;

    template <typename DirectedOption, typename VertexT, typename EdgeT, typename EdgeLengthT = double, typename VertexLabelT = std::string>
    class PhylogeneticNetwork;

    template <typename DirectedOption, typename VertexT, typename EdgeT, typename EdgeLengthT, typename VertexLabelT>
    class PhylogeneticNetworkVertex : public Detail::GraphVertex<VertexT, std::string, PhylogeneticNetwork<DirectedOption, VertexT, EdgeT, EdgeLengthT, VertexLabelT>, PhylogeneticNetworkVertex<DirectedOption, VertexT, EdgeT, EdgeLengthT, VertexLabelT>, PhylogeneticNetworkEdge<DirectedOption, VertexT, EdgeT, EdgeLengthT, VertexLabelT>, DirectedOption >
    {
        typedef Detail::GraphVertex<VertexT, std::string, PhylogeneticNetwork<DirectedOption, VertexT, EdgeT, EdgeLengthT, VertexLabelT>, PhylogeneticNetworkVertex<DirectedOption, VertexT, EdgeT, EdgeLengthT, VertexLabelT>, PhylogeneticNetworkEdge<DirectedOption, VertexT, EdgeT, EdgeLengthT, VertexLabelT>, DirectedOption> Base;
            
    public:
        
        using Base::Base;
        
        const std::string & label() const
        {
            return Base::m_customData;
        }
        
        void setLabel(const std::string & label)
        {
            Base::m_customData = label;
        }
        
        void setLabel(std::string && label)
        {
            Base::m_customData = std::move(label);
        }        
    };

    template <typename DirectedOption, typename VertexT, typename EdgeT, typename EdgeLengthT, typename VertexLabelT>
    class PhylogeneticNetworkEdge : public Detail::GraphEdge<EdgeT, EdgeLengthT, PhylogeneticNetwork<DirectedOption, VertexT, EdgeT, EdgeLengthT, VertexLabelT>, PhylogeneticNetworkVertex<DirectedOption, VertexT, EdgeT, EdgeLengthT, VertexLabelT>, PhylogeneticNetworkEdge<DirectedOption, VertexT, EdgeT, EdgeLengthT, VertexLabelT>, DirectedOption>
    {
        typedef Detail::GraphEdge<EdgeT, EdgeLengthT, PhylogeneticNetwork<DirectedOption, VertexT, EdgeT, EdgeLengthT, VertexLabelT>, PhylogeneticNetworkVertex<DirectedOption, VertexT, EdgeT, EdgeLengthT, VertexLabelT>, PhylogeneticNetworkEdge<DirectedOption, VertexT, EdgeT, EdgeLengthT, VertexLabelT>, DirectedOption> Base;
        
    public:
        
        using Base::Base;
        
        EdgeLengthT length() const
        {
            return Base::m_customData;
        }
        
        void setLength(EdgeLengthT d) 
        {
            Base::m_customData = d;
        }
    };

    template <typename DirectedOption, typename VertexT, typename EdgeT, typename EdgeLengthT, typename VertexLabelT>
    class PhylogeneticNetwork: public Detail::Graph<PhylogeneticNetwork<DirectedOption, VertexT, EdgeT, EdgeLengthT, VertexLabelT>, PhylogeneticNetworkVertex<DirectedOption, VertexT, EdgeT, EdgeLengthT, VertexLabelT>, PhylogeneticNetworkEdge<DirectedOption, VertexT, EdgeT, EdgeLengthT, VertexLabelT>, DirectedOption>
    {
        typedef Detail::Graph<PhylogeneticNetwork<DirectedOption, VertexT, EdgeT, EdgeLengthT, VertexLabelT>, PhylogeneticNetworkVertex<DirectedOption, VertexT, EdgeT, EdgeLengthT, VertexLabelT>, PhylogeneticNetworkEdge<DirectedOption, VertexT, EdgeT, EdgeLengthT, VertexLabelT>, DirectedOption> Base;
        
    public:

        PhylogeneticNetwork()
        : Base()
        {
        }
        
        PhylogeneticNetwork(const PhylogeneticNetwork<DirectedOption, VertexT, EdgeT, EdgeLengthT, VertexLabelT> & other)
        : Base(static_cast<const Base&>(other))
        {
        }
        
        PhylogeneticNetwork(PhylogeneticNetwork<DirectedOption, VertexT, EdgeT, EdgeLengthT, VertexLabelT> && other)
        : Base(std::move(other))
        {
        }
        
        PhylogeneticNetwork & operator=(const PhylogeneticNetwork & other)
        {
            if (this != &other) {
                static_cast<Base&>(*this) = static_cast<const Base&>(other);
            }
            return *this;
        }
        
        PhylogeneticNetwork & operator=(PhylogeneticNetwork && other)
        {
            if (this != &other) {
                static_cast<Base&>(*this) = std::move(other);
            }
            return *this;
        }
    };  

    template <typename Format, typename DirectedOption, typename VertexT, typename EdgeT, typename EdgeLengthT, typename VertexLabelT>
    struct GraphParserMapping<Format, PhylogeneticNetwork<DirectedOption, VertexT, EdgeT, EdgeLengthT, VertexLabelT> >
    {
        typedef PhylogeneticNetwork<DirectedOption, VertexT, EdgeT, EdgeLengthT, VertexLabelT> MyGraph;
        typedef typename MyGraph::Vertex* VertexIdentifier;
        typedef typename MyGraph::Edge* EdgeIdentifier;
        
        VertexIdentifier addVertex(MyGraph & graph, const std::vector<std::string> & comments) {return graph.addVertex();}       
        EdgeIdentifier addEdge(MyGraph & graph, VertexIdentifier from, VertexIdentifier & to) {return graph.addEdge(from, to);}      
        EdgeIdentifier trailingEdge(MyGraph & graph) {return 0;}
        void setLabel(MyGraph & g, VertexIdentifier n, std::string && label, const std::vector<std::string> & comments) {n->setLabel(std::move(label));}
        void setLength(MyGraph & g, EdgeIdentifier e, double length, const std::vector<std::string> & comments) {if (e != 0) e->setLength(length);}
        void setHybridVertexType(MyGraph & g, VertexIdentifier n, const std::string & label, const std::vector<std::string> & comments) {}
        void setRoot(MyGraph & g, VertexIdentifier n) {}
        void setUnrooted(MyGraph & g) {}
        void abort(MyGraph & g) {g.clear();}
    };  
}

#endif // GULO_PHYLOGENETICNETWORK_HPP 