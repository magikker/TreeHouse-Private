#ifndef GULO_GRAPH_VERTEX_HPP
#define GULO_GRAPH_VERTEX_HPP

#include "exception.hpp"
#include "elementdata.hpp"
#include "graphaccessor.hpp"
#include "forward.hpp"

#include <type_traits>
#include <iostream>

namespace Gulo
{
    namespace Detail
    {
        template <typename DataT, typename CustomDataT, typename DerivedGraph, typename DerivedVertex, typename DerivedEdge, typename DirectedOption>
        class GraphVertex : public GraphElementData<DataT, CustomDataT>,
                            protected GraphAccessor,
                            public VertexExtension<DerivedEdge, DirectedOption>
        {
            friend class GraphElementAccessor;
            
        public:

            typedef DerivedGraph  Graph;
            typedef DerivedEdge   Edge;
            typedef DataT         DataType;
            typedef CustomDataT   CustomDataType;

            template <typename Data = DataType, typename CustomData = CustomDataType>
            GraphVertex(typename std::enable_if<std::is_same<Data, void>::value && std::is_same<CustomData, void>::value, DerivedGraph*>::type graph, size_t id)
            : GraphElementData<DataType, CustomDataType>(),
              m_id(id),
              m_graph(graph)
            {
                if (graph == nullptr) {
                    throw BadGraphElementException("Attempt to construct vertex associated with a null graph");
                }
                increaseGraphVertexCount(m_graph);
            }

			
			template <typename Data = DataType, typename CustomData = CustomDataType>
			GraphVertex(typename std::enable_if<!std::is_same<Data, void>::value && std::is_same<CustomData, void>::value, DerivedGraph*>::type graph, size_t id, const Data & data = DataType())
            : GraphElementData<DataType, CustomDataType>(data),
              m_id(id),
              m_graph(graph)
            {
                if (graph == nullptr) {
                    throw BadGraphElementException("Attempt to construct vertex associated with a null graph");
                }
                increaseGraphVertexCount(m_graph);
            }

			template <typename Data = DataType, typename CustomData = CustomDataType>
			GraphVertex(typename std::enable_if<std::is_same<Data, void>::value && !std::is_same<CustomData, void>::value, DerivedGraph*>::type graph, size_t id, const CustomData & customData = CustomData())
            : GraphElementData<DataType, CustomDataType>(customData),
              m_id(id),
              m_graph(graph)
            {
                if (graph == nullptr) {
                    throw BadGraphElementException("Attempt to construct vertex associated with a null graph");
                }
                increaseGraphVertexCount(m_graph);
            }

			template <typename Data = DataType, typename CustomData = CustomDataType>
			GraphVertex(typename std::enable_if<!std::is_same<Data, void>::value && !std::is_same<CustomData, void>::value, DerivedGraph*>::type graph, size_t id, const Data & data = DataType(), const CustomData & customData = CustomDataType())
           : GraphElementData<DataType, CustomDataType>(data, customData),
              m_id(id),
              m_graph(graph)
            {
                if (graph == nullptr) {
                    throw BadGraphElementException("Attempt to construct vertex associated with a null graph");
                }
                increaseGraphVertexCount(m_graph);
            }

            GraphVertex(const GraphVertex & other);
        
            GraphVertex(GraphVertex && other);

            virtual ~GraphVertex() {decreaseGraphVertexCount(m_graph);}

            GraphVertex & operator =(const GraphVertex & other);
        
            GraphVertex & operator =(GraphVertex && other);

            size_t id() const            {return m_id;}

            Graph * graph()              {return m_graph;}
            const Graph * graph() const  {return m_graph;}

            template <typename Visitor>
            void visit(Visitor && visitor)          {std::remove_reference<Visitor>::type::Traversal::template visit<DirectedOption>(static_cast<DerivedVertex*>(this), visitor);}

            template <typename Visitor>
            void visit(Visitor && visitor) const    {std::remove_reference<Visitor>::type::Traversal::template visit<DirectedOption>(static_cast<const DerivedVertex*>(this), visitor);}

            void disconnect()                       {ModifyVertexStorage<DirectedOption>::disconnect(static_cast<DerivedVertex*>(this));}
            void contract(DerivedVertex * other)    {m_graph->addEdge(static_cast<DerivedVertex*>(this), other)->contract(static_cast<DerivedVertex*>(this));}
            DerivedVertex* split()                  {return GraphContraction<DirectedOption>::splitVertex(static_cast<DerivedVertex*>(this));}                 

        protected:

            size_t m_id;
            DerivedGraph * m_graph;
        };
    }
}

#endif // GULO_GRAPH_VERTEX_HPP