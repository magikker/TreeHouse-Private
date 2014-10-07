#ifndef GULO_GRAPH_EDGE_HPP
#define GULO_GRAPH_EDGE_HPP

#include "exception.hpp"
#include "elementdata.hpp"
#include "graphaccessor.hpp"
#include "forward.hpp"

#include <type_traits>

namespace Gulo
{
    namespace Detail
    {
        template <typename DataT, typename CustomDataT, typename DerivedGraph, typename DerivedVertex, typename DerivedEdge, typename DirectedOption>
        class GraphEdge : public GraphElementData<DataT, CustomDataT>,
                          protected GraphAccessor,
                          public EdgeExtension<DerivedVertex, DerivedEdge, DirectedOption>,
                          public ParallelEdgeDetection<DerivedEdge, DirectedOption>
        {
            friend class GraphElementAccessor;
            friend class EdgeExtension<DerivedVertex, DerivedEdge, DirectedOption>;
            
        public:

            typedef DerivedGraph  Graph;
            typedef DerivedVertex Vertex;
            typedef DataT         DataType;
            typedef CustomDataT   CustomDataType;

            template <typename Data = DataType, typename CustomData = CustomDataType>
			GraphEdge(typename std::enable_if<std::is_same<Data, void>::value && std::is_same<CustomData, void>::value, Vertex*>::type v1, Vertex * v2, size_t id)
            : GraphElementData<DataType, CustomDataType>(),
              m_id(id),
              m_vertex1(v1),
              m_vertex2(v2)
            {
                construct();
            }

			template <typename Data = DataType, typename CustomData = CustomDataType>
			GraphEdge(typename std::enable_if<!std::is_same<Data, void>::value && std::is_same<CustomData, void>::value, Vertex*>::type v1, Vertex * v2, size_t id, const Data & data = Data())
            : GraphElementData<DataType, CustomDataType>(data),
              m_id(id),
              m_vertex1(v1),
              m_vertex2(v2)
            {
                construct();
            }

			template <typename Data = DataType, typename CustomData = CustomDataType>
			GraphEdge(typename std::enable_if<std::is_same<Data, void>::value && !std::is_same<CustomData, void>::value, Vertex*>::type v1, Vertex * v2, size_t id, const CustomData & customData = CustomData())
            : GraphElementData<DataType, CustomDataType>(customData),
              m_id(id),
              m_vertex1(v1),
              m_vertex2(v2)
            {
                construct();
            }

			template <typename Data = DataType, typename CustomData = CustomDataType>
			GraphEdge(typename std::enable_if<!std::is_same<Data, void>::value && !std::is_same<CustomData, void>::value, Vertex*>::type v1, Vertex * v2, size_t id, const Data & data = Data(), const CustomData & customData = CustomData())
            : GraphElementData<DataType, CustomDataType>(data, customData),
              m_id(id),
              m_vertex1(v1),
              m_vertex2(v2)
            {
                construct();
            }

            GraphEdge(const GraphEdge & other);
        
            GraphEdge(GraphEdge && other);

            virtual ~GraphEdge() {decreaseGraphEdgeCount(m_graph);}

            GraphEdge & operator =(const GraphEdge & other);
        
            GraphEdge & operator =(GraphEdge && other);

            size_t id() const                                 {return m_id;}
            
            Graph * graph()                                   {return m_graph;}
            const Graph * graph() const                       {return m_graph;}
            
            Vertex * vertex1()                                {return m_vertex1;}
            const Vertex * vertex1() const                    {return m_vertex1;}
            
            Vertex * vertex2()                                {return m_vertex2;}
            const Vertex * vertex2() const                    {return m_vertex2;}
            
            Vertex * opposite(Vertex * n)                     {return (n == m_vertex1) ? m_vertex2 : ((n == m_vertex2) ? m_vertex1 : nullptr);}
            const Vertex * opposite(const Vertex * n) const   {return (n == m_vertex1) ? m_vertex2 : ((n == m_vertex2) ? m_vertex1 : nullptr);}
            
            bool isLoop() const                               {return (m_vertex1 == m_vertex2);}

            void setVertex1(DerivedVertex * v1);
            void setVertex2(DerivedVertex * v2);

            std::pair<DerivedVertex*, DerivedEdge*> split()   {return GraphContraction<DirectedOption>::template splitEdge<DerivedVertex, DerivedEdge>(static_cast<DerivedEdge*>(this));}
            void contract(DerivedVertex * retainedVertex);

        protected:

            void construct()
            {
                if (m_vertex1 == nullptr || m_vertex2 == nullptr) {
                    throw BadGraphElementException("Attempt to construct edge from one ore more null vertices");
                } else if (m_vertex1->graph() != m_vertex2->graph()) {
                    throw BadGraphElementException("Attempt to construct edge from vertices that are members of different graphs");
                }
                m_graph = m_vertex1->graph();
                increaseGraphEdgeCount(m_graph);
                ModifyVertexStorage<DirectedOption>::addEdgeToVertexStorage(static_cast<DerivedEdge*>(this));
            }

            size_t   m_id;
            Graph *  m_graph;
            Vertex * m_vertex1;
            Vertex * m_vertex2;
        };

        template <typename DataT, typename CustomDataT, typename DerivedGraph, typename DerivedVertex, typename DerivedEdge, typename DirectedOption>
        void GraphEdge<DataT, CustomDataT, DerivedGraph, DerivedVertex, DerivedEdge, DirectedOption>::setVertex1(DerivedVertex * v1)
        {
            if (v1 == nullptr || v1->graph() != this->graph()) {
                throw Gulo::BadGraphElementException("Attempt to set vertex 1 of an edge to a vertex that is not a member of the edge's graph");
            }
            if (this->m_vertex1 != v1) {
                ModifyVertexStorage<DirectedOption>::setVertex1Storage(static_cast<DerivedEdge*>(this), v1);
                m_vertex1 = v1;    
            }
        }

        template <typename DataT, typename CustomDataT, typename DerivedGraph, typename DerivedVertex, typename DerivedEdge, typename DirectedOption>
        void GraphEdge<DataT, CustomDataT, DerivedGraph, DerivedVertex, DerivedEdge, DirectedOption>::setVertex2(DerivedVertex * v2)
        {
            if (v2 == nullptr || v2->graph() != this->graph()) {
                throw Gulo::BadGraphElementException("Attempt to set vertex 2 of an edge to a vertex that is not a member of the edge's graph");
            }
            if (this->m_vertex2 != v2) {
                ModifyVertexStorage<DirectedOption>::setVertex2Storage(static_cast<DerivedEdge*>(this), v2); 
                m_vertex2 = v2;
            }
        }

        template <typename DataT, typename CustomDataT, typename DerivedGraph, typename DerivedVertex, typename DerivedEdge, typename DirectedOption>
        void GraphEdge<DataT, CustomDataT, DerivedGraph, DerivedVertex, DerivedEdge, DirectedOption>::contract(DerivedVertex * retainedVertex = nullptr)
        {
            if (retainedVertex != nullptr && retainedVertex != m_vertex1 && retainedVertex != m_vertex2) {
                throw BadGraphElementException("Attempt to contract an edge, retaining a vertex that is not a member of the edge");
            }
            if (isLoop()) {
                m_graph->removeEdge(static_cast<DerivedEdge*>(this));
            } else {
                if (retainedVertex == 0 || retainedVertex == m_vertex1) {
                    GraphContraction<DirectedOption>::contractEdge(static_cast<DerivedEdge*>(this), m_vertex1, m_vertex2);
                } else {
                    GraphContraction<DirectedOption>::contractEdge(static_cast<DerivedEdge*>(this), m_vertex2, m_vertex1);
                }
            }
        }
    }
}

#endif // GULO_GRAPH_EDGE_HPP