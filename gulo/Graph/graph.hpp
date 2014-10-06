#ifndef GULO_GRAPH_GRAPH_HPP
#define GULO_GRAPH_GRAPH_HPP

#include <random>
#include <boost/pool/pool.hpp>
#include <algorithm>

#include "exception.hpp"
#include "elementdata.hpp"
#include "elementvector.hpp"
#include "elementvectoraccessor.hpp"
#include "elementaccessor.hpp"
#include "forward.hpp"
#include "graphparser.hpp"

namespace Gulo
{
    namespace Detail
    {      
        template <typename Vertex, typename VertexData, typename VertexCustomData>
        struct GraphAddVertex : virtual protected GraphElementVectorAccessor
        {
            Vertex * addVertex(const VertexData & data = VertexData(), const VertexCustomData & customData = VertexCustomData())
            {
                Vertex * n = new(static_cast<typename Vertex::Graph*>(this)->m_vertexPool->malloc()) Vertex(static_cast<typename Vertex::Graph*>(this), static_cast<typename Vertex::Graph*>(this)->m_vertices.size(), data, customData);
                this->getVector(static_cast<typename Vertex::Graph*>(this)->m_vertices).push_back(n);
                return n;
            }
        };

        template <typename Vertex, typename VertexData>
        struct GraphAddVertex<Vertex, VertexData, void> : virtual protected GraphElementVectorAccessor
        {
            Vertex * addVertex(const VertexData & data = VertexData())
            {
                Vertex * n = new(static_cast<typename Vertex::Graph*>(this)->m_vertexPool->malloc()) Vertex(static_cast<typename Vertex::Graph*>(this), static_cast<typename Vertex::Graph*>(this)->m_vertices.size(), data);
                this->getVector(static_cast<typename Vertex::Graph*>(this)->m_vertices).push_back(n);
                return n;
            }
        };

        template <typename Vertex, typename VertexCustomData>
        struct GraphAddVertex<Vertex, void, VertexCustomData> : virtual protected GraphElementVectorAccessor
        {
            Vertex * addVertex(const VertexCustomData & customData = VertexCustomData())
            {
                Vertex * n = new(static_cast<typename Vertex::Graph*>(this)->m_vertexPool->malloc()) Vertex(static_cast<typename Vertex::Graph*>(this), static_cast<typename Vertex::Graph*>(this)->m_vertices.size(), customData);
                this->getVector(static_cast<typename Vertex::Graph*>(this)->m_vertices).push_back(n);
                return n;
            }
        };

        template <typename Vertex>
        struct GraphAddVertex<Vertex, void, void> : virtual protected GraphElementVectorAccessor
        {
            Vertex * addVertex()
            {
                Vertex * n = new(static_cast<typename Vertex::Graph*>(this)->m_vertexPool->malloc()) Vertex(static_cast<typename Vertex::Graph*>(this), static_cast<typename Vertex::Graph*>(this)->m_vertices.size());
                this->getVector(static_cast<typename Vertex::Graph*>(this)->m_vertices).push_back(n);
                return n;
            }
        };

        template <typename Edge, typename Vertex, typename EdgeData, typename EdgeCustomData>
        struct GraphAddEdge : virtual protected GraphElementVectorAccessor
        {
            Edge * addEdge(Vertex * v1, Vertex * v2, const EdgeData & data = EdgeData(), const EdgeCustomData & customData = EdgeCustomData())
            {
                Edge * e = new(static_cast<typename Edge::Graph*>(this)->m_edgePool->malloc()) Edge(v1, v2, static_cast<typename Edge::Graph*>(this)->m_edges.size(), data, customData);
                this->getVector(static_cast<typename Edge::Graph*>(this)->m_edges).push_back(e);
                return e;
            }
        };

        template <typename Edge, typename Vertex, typename EdgeData>
        struct GraphAddEdge<Edge, Vertex, EdgeData, void> : virtual protected GraphElementVectorAccessor
        {
            Edge * addEdge(Vertex * v1, Vertex * v2, const EdgeData & data = EdgeData())
            {
                Edge * e = new(static_cast<typename Edge::Graph*>(this)->m_edgePool->malloc()) Edge(v1, v2, static_cast<typename Edge::Graph*>(this)->m_edges.size(), data);
                this->getVector(static_cast<typename Edge::Graph*>(this)->m_edges).push_back(e);
                return e;
            }
        };

        template <typename Edge, typename Vertex,  typename EdgeCustomData>
        struct GraphAddEdge<Edge, Vertex, void, EdgeCustomData> : virtual protected GraphElementVectorAccessor
        {
            Edge * addEdge(Vertex * v1, Vertex * v2, const EdgeCustomData & customData = EdgeCustomData())
            {
                Edge * e = new(static_cast<typename Edge::Graph*>(this)->m_edgePool->malloc()) Edge(v1, v2, static_cast<typename Edge::Graph*>(this)->m_edges.size(), customData);
                this->getVector(static_cast<typename Edge::Graph*>(this)->m_edges).push_back(e);
                return e;
            }
        };

		template <typename Edge, typename Vertex>
		struct GraphAddEdge<Edge, Vertex, void, void> : virtual protected GraphElementVectorAccessor
		{
			Edge * addEdge(Vertex * v1, Vertex * v2)
			{
				Edge * e = new(static_cast<typename Edge::Graph*>(this)->m_edgePool->malloc()) Edge(v1, v2, static_cast<typename Edge::Graph*>(this)->m_edges.size());
				this->getVector(static_cast<typename Edge::Graph*>(this)->m_edges).push_back(e);
				return e;
			}
		};

        template <typename DerivedGraph, typename DerivedVertex, typename DerivedEdge, typename DirectedOption>
        class Graph : public GraphAddVertex<DerivedVertex, typename DerivedVertex::DataType, typename DerivedVertex::CustomDataType>,
                      public GraphAddEdge<DerivedEdge, DerivedVertex, typename DerivedEdge::DataType, typename DerivedEdge::CustomDataType>,
                      protected virtual GraphElementVectorAccessor,
                      protected virtual GraphElementAccessor,
                      public GraphExtension<DerivedGraph, DerivedVertex, DirectedOption>
        {      
            friend class GraphAccessor;
            typedef GraphElementData<typename DerivedVertex::DataType, typename DerivedVertex::CustomDataType> VertexElementData;
            typedef GraphElementData<typename DerivedEdge::DataType, typename DerivedEdge::CustomDataType> EdgeElementData;
            template <typename, typename, typename> friend struct GraphAddVertex;
            template <typename, typename, typename, typename> friend struct GraphAddEdge;

        public:

            typedef DerivedVertex                           Vertex;
            typedef DerivedEdge                             Edge;
            typedef typename DerivedVertex::DataType        VertexDataType;
            typedef typename DerivedVertex::CustomDataType  VertexCustomDataType;
            typedef typename DerivedEdge::DataType          EdgeDataType;
            typedef typename DerivedEdge::CustomDataType    EdgeCustomDataType;
            typedef DirectedOption                          Directionality;

            typedef boost::pool<boost::default_user_allocator_malloc_free> MemoryPool;
            
            Graph() 
            : m_pausedIndexing(false),
              m_vertexCount(0),
              m_edgeCount(0),
              m_vertexPool(new MemoryPool(sizeof(DerivedVertex))),
              m_edgePool(new MemoryPool(sizeof(DerivedEdge)))
            {
            } 

            Graph(const Graph<DerivedGraph, DerivedVertex, DerivedEdge, DirectedOption> & other)
            : Graph()
            {
                *this = other;
            }
             
            Graph(Graph<DerivedGraph, DerivedVertex, DerivedEdge, DirectedOption> && other)
            : m_pausedIndexing(other.m_pausedIndexing),
              m_vertexCount(other.m_vertexCount),
              m_edgeCount(other.m_edgeCount),
              m_vertices(std::move(other.m_vertices)),
              m_edges(std::move(other.m_edges)),
              m_vertexPool(other.m_vertexPool),
              m_edgePool(other.m_edgePool)
            {
                other.m_pausedIndexing = false;
                other.m_vertexCount = 0;
                other.m_edgeCount = 0;
                other.m_vertexPool = nullptr;
                other.m_edgePool = nullptr;
            }
             
            virtual ~Graph() 
            {
                clear();
                if (m_vertexPool != nullptr) {
                    delete m_vertexPool;
                }
                if (m_edgePool != nullptr) {
                    delete m_edgePool;
                }
            }

            Graph & operator=(const Graph & other)
            {
                if (this != &other) {
                    clear();
                    m_pausedIndexing = other.m_pausedIndexing;
                    m_vertices.reserve(other.m_vertices.size());
                    m_edges.reserve(other.m_edges.size());
                    for (auto * v : other.m_vertices) {
                        if (v != nullptr) {
                            getVector(m_vertices).push_back(new(m_vertexPool->malloc()) DerivedVertex(static_cast<typename DerivedVertex::Graph*>(this), v->id()));
                            static_cast<VertexElementData&>(*(m_vertices.back())) = static_cast<const VertexElementData&>(*v);
                            CopyVertexExtension<DirectedOption>::template copy<Vertex>(m_vertices.back(), v);
                        } else {
                            getVector(m_vertices).push_back(0);
                        }
                    }                    
                    for (auto * e : other.m_edges) {
                        if (e != nullptr) {
                            auto * edge = (new(m_edgePool->malloc()) Edge(vertex(e->vertex1()->id()), vertex(e->vertex2()->id()), e->id()));
                            static_cast<EdgeElementData&>(*edge) = static_cast<const EdgeElementData&>(*e);
                            CopyEdgeExtension<DirectedOption>::template copy<Edge>(edge, e);
                            getVector(m_edges).push_back(edge);
                        } else {
                            getVector(m_edges).push_back(nullptr);
                        }
                    }
                }
                return *this;
            }
             
            Graph & operator=(Graph && other)
            {
                clear();
                if (m_vertexPool != nullptr) {
                    delete m_vertexPool;
                }
                if (m_edgePool != nullptr) {
                    delete m_edgePool;
                }
                m_pausedIndexing = other.m_pausedIndexing;  
                m_vertexCount = other.m_vertexCount; 
                m_edgeCount = other.m_edgeCount;  
                m_vertices = std::move(other.m_vertices); 
                m_edges = std::move(other.m_edges);  
                m_vertexPool = other.m_vertexPool;              
                m_edgePool = other.m_edgePool;
                other.m_pausedIndexing = false;
                other.m_vertexCount = 0;
                other.m_edgeCount = 0;
                other.m_vertexPool = nullptr;
                other.m_edgePool = nullptr;
                return *this;
            }

            template <typename Format, typename Mapping = GraphParserMapping<Format, DerivedGraph>, typename Data>
            static DerivedGraph parse(const Data & data, typename Format::template Parser<Mapping>::ParseResult * result = 0)
            {
                DerivedGraph g;
                GraphParser<Format, DerivedGraph, Mapping> parser;
                if (result != 0) {
                    *result = parser.parse(data, g);
                } else {
                    parser.parse(data, g);
                }
                return g;
            }

            DerivedGraph subgraph(const Vertex * vertex) const
            {
                struct Visitor : public DepthFirstVisitor
                {
                    Visitor(size_t vertexCount)
                    {
                        activeEdge = nullptr;
                        newVertices.assign(vertexCount, nullptr);
                    }

                    void reachUnvisitedVertex(const Vertex * v)
                    {
                        newVertices[v->id()] = graph.addVertex();
                        static_cast<VertexElementData&>(*(newVertices[v->id()])) = static_cast<const VertexElementData&>(*v);
                        CopyVertexExtension<DirectedOption>::template copy<Vertex>(newVertices[v->id()], v);
                        if (activeEdge != nullptr) {
                            auto * edge = graph.addEdge(newVertices[activeEdge->vertex1()->id()], newVertices[activeEdge->vertex2()->id()]);
                            static_cast<EdgeElementData&>(*edge) = static_cast<const EdgeElementData&>(*activeEdge);
                            CopyEdgeExtension<DirectedOption>::template copy<Edge>(edge, activeEdge);
                        }
                    }

                    void startEdge(const Edge * e)
                    {
                        activeEdge = e;
                    }

                    DerivedGraph graph;
                    std::vector<Vertex*> newVertices;
                    const Edge * activeEdge;
                };
                Visitor vis(m_vertices.size());
                vertex->visit(vis);
                return vis.graph;
            }
            
            size_t vertexCount() const                              {return m_vertexCount;}
            size_t edgeCount() const                                {return m_edgeCount;} 

            Vertex * vertex(size_t i)                               {return m_vertices[i];} 
            const Vertex * vertex(size_t i) const                   {return m_vertices[i];}
            Edge * edge(size_t i)                                   {return m_edges[i];}
            const Edge * edge(size_t i) const                       {return m_edges[i];}
             
            GraphElementVector<Vertex> & vertices()                 {return m_vertices;}
             
            const GraphElementVector<Vertex> & vertices() const     {return m_vertices;}
            GraphElementVector<DerivedEdge> & edges()               {return m_edges;}
            const GraphElementVector<DerivedEdge> & edges() const   {return m_edges;}

            template <typename RNG>
            Vertex * randomVertex(RNG & rng)                        
            {
                std::uniform_int_distribution<size_t> gen(0, m_vertices.size()-1);
                return m_vertices[gen(rng)];
            }
             
            template <typename RNG>
            const DerivedVertex * randomVertex(RNG & rng) const
            {
                std::uniform_int_distribution<size_t> gen(0, m_vertices.size()-1);
                return m_vertices[gen(rng)];
            }
            
            template <typename RNG>
            DerivedEdge * randomEdge(RNG & rng)
            {
                std::uniform_int_distribution<size_t> gen(0, m_edges.size()-1);
                return m_edges[gen(rng)];
            }
             
            template <typename RNG>
            const DerivedEdge * randomEdge(RNG & rng) const
            {
                std::uniform_int_distribution<size_t> gen(0, m_edges.size()-1);
                return m_edges[gen(rng)];
            }

            void clear()
            {
                for (auto * v : m_vertices) {
                    if (v != nullptr) {
                        v->~DerivedVertex();
                        m_vertexPool->free(v);
                    }
                }
                for (auto * e : m_edges) {
                    if (e != nullptr) {
                        e->~DerivedEdge();
                        m_edgePool->free(e);
                    }
                }
                m_pausedIndexing = false;
                m_vertexCount = 0;
                m_edgeCount = 0;
                getVector(m_vertices).clear();
                getVector(m_edges).clear(); 
            }

            void removeVertex(DerivedVertex * v)
            {
                if (v == nullptr) {
                    throw BadGraphElementException("Attempt to remove a null vertex");
                }
                if (v->graph() != this) {
                    throw BadGraphElementException("Attempt to remove a vertex that is not a member of the graph");
                }
                v->disconnect();
                auto & vertexVec = this->getVector(m_vertices);
                if (!elementReindexingIsSuspended()) {
                    std::for_each(vertexVec.begin() + v->id() + 1, vertexVec.end(), [&](DerivedVertex * vv)->void{--(getElementId(vv));});
                    vertexVec.erase(vertexVec.begin() + v->id());
                } else {
                    *(vertexVec.begin() + v->id()) = nullptr;
                }
                v->~DerivedVertex();
                m_vertexPool->free(v);
            }

            void removeEdge(DerivedEdge * e)
            {
                if (e == nullptr) {
                    throw BadGraphElementException("Attempt to remove a null edge");
                }
                if (e->graph() != this) {
                    throw BadGraphElementException("Attempt to remove an edge that is not a member of the graph");
                }
                ModifyVertexStorage<DirectedOption>::removeEdgeFromVertexStorage(e);
                auto & edgeVec = getVector(m_edges);
                if (!elementReindexingIsSuspended()) {
                    std::for_each(edgeVec.begin() + e->id() + 1, edgeVec.end(), [&](DerivedEdge * e)->void{--(getElementId(e));});
                    edgeVec.erase(edgeVec.begin() + e->id());
                } else {
                    *(edgeVec.begin() + e->id()) = nullptr;
                }
                e->~DerivedEdge();
                this->m_edgePool->free(e);
            }

            void suspendElementReindexing()
            {
                m_pausedIndexing = true;
            }
             
            bool elementReindexingIsSuspended() const
            {
                return m_pausedIndexing;
            }
             
            void resumeElementReindexing()
            {
                if (elementReindexingIsSuspended()) {
                    m_pausedIndexing = false;
                    auto & edgeVec = getVector(m_edges);
                    auto & vertVec = getVector(m_vertices);
                    edgeVec.erase(std::remove(edgeVec.begin(), edgeVec.end(), nullptr), edgeVec.end());
                    vertVec.erase(std::remove(vertVec.begin(), vertVec.end(), nullptr), vertVec.end());
                    for (size_t i = 0; i < m_edges.size(); ++i) { 
                        getElementId(m_edges[i]) = i;
                    }
                    for (size_t i = 0; i < m_vertices.size(); ++i) {
                        getElementId(m_vertices[i]) = i;
                    }
                }
            }

            template <typename Visitor>
            void visit(Visitor && visitor)
            {
                GraphTraversal<DirectedOption>::template visit<DerivedGraph>(static_cast<DerivedGraph&>(*this), visitor);
            }

            template <typename Visitor>
            void visit(Visitor && visitor) const
            {
                GraphTraversal<DirectedOption>::template visit<DerivedGraph>(static_cast<const DerivedGraph&>(*this), visitor);
            }

            void simplify() // todo, some options to allow merging of parallel edges or loops prior to removal?
            {
                bool b = elementReindexingIsSuspended();
                if (!b) {
                    suspendElementReindexing();
                }
                std::for_each(m_edges.begin(), m_edges.end(), [&](Edge * e)->void{if (e->isLoop() || e->existsParallelEdge()) ModifyVertexStorage<DirectedOption>::removeEdgeFromVertexStorage(e);});
                if (!b) {
                    resumeElementReindexing();
                }
            }

            DerivedGraph disjunctUnion(const DerivedGraph & other) const
            {
                DerivedGraph g = static_cast<const DerivedGraph&>(*this);
                g.disjunctUniteWith(other);
                return g;
            }

            DerivedGraph join(const DerivedGraph & other, const std::function<void(DerivedEdge*)> & setEdgeData = [](DerivedEdge*)->void{}) const
            {
                DerivedGraph g = static_cast<const DerivedGraph&>(*this);
                g.template joinWith(other, setEdgeData);
                return g;
            }

            DerivedGraph sum(const DerivedGraph & other, const std::function<bool(const DerivedVertex * v1, const DerivedVertex * v2, const DerivedVertex * w1, const DerivedVertex * w2, const DerivedEdge * e, const DerivedEdge * f)> & test = StandardEdgeIsomorphismTest<DirectedOption>()) const
            {
                DerivedGraph g = static_cast<const DerivedGraph&>(*this);
                g.template sumWith(other, test);
                return g;
            }

            DerivedGraph difference(const DerivedGraph & other, const std::function<bool(const DerivedVertex * v1, const DerivedVertex * v2, const DerivedVertex * w1, const DerivedVertex * w2, const DerivedEdge * e, const DerivedEdge * f)> & test = StandardEdgeIsomorphismTest<DirectedOption>()) const
            {
                DerivedGraph g = static_cast<const DerivedGraph&>(*this);
                g.template differenceWith(other, test);
                return g;
            }

            DerivedGraph intersection(const DerivedGraph & other, const std::function<bool(const DerivedVertex * v1, const DerivedVertex * v2, const DerivedVertex * w1, const DerivedVertex * w2, const DerivedEdge * e, const DerivedEdge * f)> & test = StandardEdgeIsomorphismTest<DirectedOption>()) const
            {
                DerivedGraph g = static_cast<const DerivedGraph&>(*this);
                g.template intersectWith(other, test);
                return g;
            }

            DerivedGraph cartesianProduct      (const DerivedGraph & other, const std::function<void(DerivedVertex * target, const DerivedVertex * v1, const DerivedVertex * v2)> & vertexMultiplier = [](DerivedVertex*, const DerivedVertex*, const DerivedVertex*)->void{}) const;
            DerivedGraph rootedProduct         (const DerivedGraph & other, const Vertex * root, const std::function<void(DerivedVertex * target, const DerivedVertex * v1, const DerivedVertex * v2)> & vertexMultiplier = [](DerivedVertex*, const DerivedVertex*, const DerivedVertex*)->void{}) const;
            DerivedGraph lexicographicalProduct(const DerivedGraph & other, const std::function<void(DerivedVertex * target, const DerivedVertex * v1, const DerivedVertex * v2)> & vertexMultiplier = [](DerivedVertex*, const DerivedVertex*, const DerivedVertex*)->void{}) const;
            DerivedGraph strongProduct         (const DerivedGraph & other, const std::function<void(DerivedVertex * target, const DerivedVertex * v1, const DerivedVertex * v2)> & vertexMultiplier = [](DerivedVertex*, const DerivedVertex*, const DerivedVertex*)->void{}, const std::function<void(DerivedEdge * target, const DerivedEdge * e1, const DerivedEdge * e2)> & edgeMultiplier = [](DerivedEdge*, const DerivedEdge*, const DerivedEdge*)->void{}) const;
            DerivedGraph tensorProduct         (const DerivedGraph & other, const std::function<void(DerivedVertex * target, const DerivedVertex * v1, const DerivedVertex * v2)> & vertexMultiplier = [](DerivedVertex*, const DerivedVertex*, const DerivedVertex*)->void{}, const std::function<void(DerivedEdge * target, const DerivedEdge * e1, const DerivedEdge * e2)> & edgeMultiplier = [](DerivedEdge*, const DerivedEdge*, const DerivedEdge*)->void{}) const;


            bool isAcyclicGraph() const
            {
                struct GraphBaseAcyclicTester : public DepthFirstVisitor
                {
                    void reachQueuedVertex(const Vertex * v) {isAcyclic = false;}
                    bool terminate() const {return !isAcyclic;}
                    bool isAcyclic{true};
                } t;
                this->visit(t);
                return t.isAcyclic;
            }

        protected:

            void joinWith(const DerivedGraph & other, 
                          const std::function<void(DerivedEdge*)> & setEdgeData = [](DerivedEdge*)->void{})
            {
                size_t vc = m_vertices.size();
                disjunctUniteWith(other);
                size_t vc2 = m_vertices.size();
                for (size_t i = 0; i < vc; ++i) {
                    for (size_t j = vc; j < vc2; ++j) {
                        setEdgeData(this->addEdge(vertex(i), vertex(j)));
                    }
                }
            }           

            void disjunctUniteWith(const DerivedGraph & other)
            {
                size_t vertexOffset = m_vertices.size();
                DerivedGraph * graph = static_cast<DerivedGraph*>(this);
                m_vertices.reserve(m_vertices.size() + other.m_vertices.size());
                for (auto * v : other.m_vertices) {
                    if (v != 0) {
                        getVector(m_vertices).push_back(new(m_vertexPool->malloc()) DerivedVertex(graph, v->id() + vertexOffset));
                        static_cast<VertexElementData&>(*(m_vertices.back())) = static_cast<const VertexElementData&>(*v); 
                        CopyVertexExtension<DirectedOption>::template copy<Vertex>(m_vertices.back(), v);
                    } else {
                        getVector(m_vertices).push_back(0);
                    }
                }
                size_t edgeOffset = m_edges.size();
                size_t vc = other.m_edges.size();
                m_edges.reserve(m_edges.size() + other.m_edges.size());
                for (size_t i = 0; i < vc; ++i) {
                    auto * e = other.m_edges[i];
                    if (e != 0) {
                        auto * edge = (new(m_edgePool->malloc()) Edge(vertex(e->vertex1()->id() + vertexOffset), vertex(e->vertex2()->id() + vertexOffset), e->id() + edgeOffset));
                        static_cast<EdgeElementData&>(*edge) = static_cast<const EdgeElementData&>(*e);
                        CopyEdgeExtension<DirectedOption>::template copy<Edge>(edge, e);
                        getVector(m_edges).push_back(edge);
                    } else {
                        getVector(m_edges).push_back(0);
                    }
                }
            }

            void sumWith(const DerivedGraph & other, 
                         const std::function<bool(const DerivedVertex * v1, const DerivedVertex * v2, const DerivedVertex * w1, const DerivedVertex * w2, const DerivedEdge * e, const DerivedEdge * f)> & test = StandardEdgeIsomorphismTest<DirectedOption>())
            {
                for (auto * e : other.edges()) {
                    if (e != nullptr) {
                        bool add = true;
                        for (auto * k : m_vertices[e->vertex1()->id()]->edgesWith(m_vertices[e->vertex2()->id()])) {
                            if (test(m_vertices[e->vertex1()->id()], m_vertices[e->vertex2()->id()], e->vertex1(), e->vertex2(), k, e)) {
                                add = false;
                                break;
                            }
                        }
                        if (add) {
                            auto * edge = this->addEdge(m_vertices[e->vertex1()->id()], m_vertices[e->vertex2()->id()]);
                            CopyEdgeExtension<DirectedOption>::template copy<Edge>(edge, e);
                        }
                    }
                }
            }

            void intersectWith(const DerivedGraph & other, const std::function<bool(const DerivedVertex * v1, const DerivedVertex * v2, const DerivedVertex * w1, const DerivedVertex * w2, const DerivedEdge * e, const DerivedEdge * f)> & test = StandardEdgeIsomorphismTest<DirectedOption>())
            {
                bool b = elementReindexingIsSuspended();
                if (!b) {
                    suspendElementReindexing();
                }
                for (auto * e : m_edges) {
                    if (e != nullptr) {
                        bool remove = true;
                        for (auto * k : other.vertex(e->vertex1()->id())->edgesWith(other.vertex(e->vertex2()->id()))) {
                            if (test(e->vertex1(), e->vertex2(), other.vertex(e->vertex1()->id()), other.vertex(e->vertex2()->id()), e, k)) {
                                remove = false;
                                break;
                            }
                        }
                        if (remove) {
                            this->removeEdge(e);
                        }
                    }
                }
                if (!b) {
                    resumeElementReindexing();
                }
            }

            void differenceWith(const DerivedGraph & other, const std::function<bool(const DerivedVertex * v1, const DerivedVertex * v2, const DerivedVertex * w1, const DerivedVertex * w2, const DerivedEdge * e, const DerivedEdge * f)> & test = StandardEdgeIsomorphismTest<DirectedOption>())
            {
                bool b = elementReindexingIsSuspended();
                if (!b) {
                    suspendElementReindexing();
                }
                for (auto * e : m_edges) {
                    if (e != nullptr) {
                        for (auto * k : other.vertex(e->vertex1()->id())->edgesWith(other.vertex(e->vertex2()->id()))) {
                            if (test(e->vertex1(), e->vertex2(), other.vertex(e->vertex1()->id()), other.vertex(e->vertex2()->id()), e, k)) {
                                this->removeEdge(e);
                                break;
                            }
                        }
                    }
                }
                if (!b) {
                    resumeElementReindexing();
                }
            }
 
            bool m_pausedIndexing;
            size_t m_vertexCount;
            size_t m_edgeCount;
            GraphElementVector<DerivedVertex> m_vertices;
            GraphElementVector<DerivedEdge> m_edges;
            MemoryPool * m_vertexPool;
            MemoryPool * m_edgePool;
        };   
    }    
}    

#include "graphproduct.Hpp"

#endif // GULO_GRAPH_HPP