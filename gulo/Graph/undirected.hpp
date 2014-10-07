#ifndef GULO_GRAPH_EDGEUNDIRECTED_HPP
#define GULO_GRAPH_EDGEUNDIRECTED_HPP

#include "elementvectoraccessor.hpp"
#include "elementvector.hpp"
#include "forward.hpp"
#include "graphaccessor.hpp"
#include "elementaccessor.hpp"
#include "elementdata.hpp"
#include "visitorstate.hpp"
#include "graphtraversal.hpp"
#include "weaklyconnectedcomponentsundirected.hpp"

#include <boost/range/adaptor/reversed.hpp>

#include <algorithm>

namespace Gulo
{
    struct Undirected;

    namespace Detail
    {
        /////////////////////
        // Graph Extension //
        /////////////////////

        template <typename Graph, typename Vertex>
        class GraphExtension<Graph, Vertex, Undirected>
        {
        public:

            constexpr static bool isUndirectedGraph() {return true;}
            constexpr static bool isDirectedGraph()   {return false;}
            constexpr static bool isBidirectedGraph() {return false;}
            constexpr static bool storesOutEdges()    {return false;}
            constexpr static bool storesInEdges()     {return false;}  

            std::vector<std::set<Vertex*>> connectedComponents()
            {
                UndirectedWeaklyConnectedComponentGatherer<Vertex*> v;
                static_cast<Graph*>(this)->visit(v);
                return v.components;
            }

            std::vector<std::set<const Vertex*>> connectedComponents() const
            {
                UndirectedWeaklyConnectedComponentGatherer<const Vertex*> v;
                static_cast<const Graph*>(this)->visit(v);
                return v.components;
            }

            bool isConnectedGraph() const
            {
                UndirectedWeaklyConnectedGraphTester<const Vertex*> v;
                static_cast<const Graph*>(this)->visit(v);
                return !(v.terminate());
            }
        };

        //////////////////////
        // Vertex Extension //
        //////////////////////

        template <typename Edge>
        class VertexExtension<Edge, Undirected>
        {
            typedef typename Edge::Vertex Vertex;

        public:

            GraphElementVector<Edge> & edges()             {return m_edges;}
            const GraphElementVector<Edge> & edges() const {return m_edges;}

            std::vector<const Edge *> edgesWith(const Vertex * v) const
            {
                std::vector<const Edge *> result;
                if (v != 0) {
                    std::for_each(this->m_edges.begin(), this->m_edges.end(), [&](const Edge*e)->void{ if (v == e->opposite(static_cast<const Vertex*>(this))) result.push_back(e);});
                }
                return result;
            }
            
            std::vector<Edge *> edgesWith(const Vertex * v)
            {
                std::vector<Edge *> result;
                if (v != 0) {
                    std::for_each(this->m_edges.begin(), this->m_edges.end(), [&](const Edge*e)->void{ if (v == e->opposite(static_cast<const Vertex*>(this))) result.push_back(const_cast<Edge*>(e));});
                }
                return result;
            }

            bool isAdjacentTo(const Vertex * v) const
            {
                for (auto * e : this->m_edges) {
                    if (v == e->opposite(static_cast<const Vertex*>(this))) {
                       return true;
                    }
                }
                return false;
            }

            size_t degree() const
            {
                size_t result = 0;
                for (auto * e : this->m_edges) {
                    result += (e->isLoop()) ? 2 : 1;
                }
                return result;
            }

        private:

            GraphElementVector<Edge> m_edges;
        };

        /////////////////////////
        // Copy Edge Extension //
        /////////////////////////

        template <>
        struct CopyEdgeExtension<Undirected>
        {
            template <typename Edge>
            static void copy(Edge * target, const Edge * source)
            {
            }
        };

        ///////////////////////////
        // Copy Vertex Extension //
        ///////////////////////////

        template <>
        struct CopyVertexExtension<Undirected>
        {
            template <typename Edge>
            static void copy(Edge * target, const Edge * source)
            {
            }
        };

        ///////////////////////////
        // Modify Vertex Storage //
        ///////////////////////////

        template <>
        struct ModifyVertexStorage<Undirected> : protected GraphElementVectorAccessor,
                                                 protected GraphAccessor,
                                                 protected GraphElementAccessor
        {
            template <typename Edge>
            static void addEdgeToVertexStorage(Edge * e) 
            {
                getVector(e->vertex1()->edges()).push_back(e);
                if (!e->isLoop()) {
                    getVector(e->vertex2()->edges()).push_back(e);
                }
            }

            template <typename Edge>
            static void removeEdgeFromVertexStorage(Edge * e) 
            {
                if (e->vertex1() != nullptr) {
                    auto & vec = getVector(e->vertex1()->edges());
                    vec.erase(std::find(vec.begin(), vec.end(), e));
                }
                if (e->vertex2() != nullptr && e->vertex2() != e->vertex1()) {
                    auto & vec = getVector(e->vertex2()->edges());
                    vec.erase(std::find(vec.begin(), vec.end(), e));
                }
            }

            template <typename Edge, typename Vertex>
            static void setVertex1Storage(Edge * e, Vertex * v)
            {
                if (!e->isLoop() || e->vertex2() == v) {
                    auto & vec = getVector(e->vertex1()->edges());
                    vec.erase(std::find(vec.begin(), vec.end(), e));
                }
                if (e->vertex2() != v) {
                    getVector(v->edges()).push_back(e);
                }
            }

            template <typename Edge, typename Vertex>
            static void setVertex2Storage(Edge * e, Vertex * v)
            {
                if (!e->isLoop() || e->vertex1() == v) {
                    auto & vec = getVector(e->vertex2()->edges());
                    vec.erase(std::find(vec.begin(), vec.end(), e));
                }
                if (e->vertex1() != v) {
                    getVector(v->edges()).push_back(e);
                }
            }

            template <typename Vertex>
            static void disconnect(Vertex * v)
            {
                typedef typename Vertex::Edge Edge;
                auto & edgeVec = getVector(v->graph()->edges());
                size_t minEdgeIndex = v->graph()->edges().size();
                for (auto * e : v->edges()) {
                    minEdgeIndex = std::min(e->id(), minEdgeIndex);
                    auto & E = getVector(e->opposite(v)->edges());
                    E.erase(std::find(E.begin(), E.end(), e));
                    if (!e->isLoop()) {
                        edgeVec[e->id()] = 0;
                        e->~Edge();
                        edgePool(v->graph())->free(e);
                    } else {
                        if (edgeVec[e->id()] == 0) {
                            e->~Edge();
                            edgePool(v->graph())->free(e);
                        } else {
                            edgeVec[e->id()] = 0;
                        }
                    }
                }
                if (!v->graph()->elementReindexingIsSuspended()) {
                    for (size_t i = minEdgeIndex; i < edgeVec.size(); ++i) {
                        if (edgeVec[i] == 0) {
                            edgeVec.erase(edgeVec.begin() + i);
                            --i;
                        } else {
                            getElementId(edgeVec[i]) = i;
                        }
                    }
                }
                getVector(v->edges()).clear();
            }
        };

        /////////////////////////////
        // Parallel Edge Detection //
        /////////////////////////////

        template <typename Edge>
        struct ParallelEdgeDetection<Edge, Undirected>
        {
            bool existsParallelEdge() const
            {
                for (auto * e : (static_cast<const Edge*>(this)->vertex1()->edges().size() < static_cast<const Edge*>(this)->vertex2()->edges().size()) ? static_cast<const Edge*>(this)->vertex1()->edges() : static_cast<const Edge*>(this)->vertex2()->edges()) {
                    if (e != this && ((e->vertex1() == static_cast<const Edge*>(this)->vertex1() && e->vertex2() == static_cast<const Edge*>(this)->vertex2()) || (e->vertex1() == static_cast<const Edge*>(this)->vertex2() && e->vertex2() == static_cast<const Edge*>(this)->vertex1()))) {
                        return true;
                    }
                }
                return false;
            }
        };

        /////////////////////////////
        // Traverse Outgoing Edges //
        /////////////////////////////

        template <>
        struct TraverseOutgoingEdges<Undirected>
        {
            template <typename Vertex, typename Func>
            static void go(Vertex* v, Func & func, const typename Vertex::Graph::Edge * prev)
            {
                for (auto * edge : v->edges()) {
                    if (prev != edge) {
                        func(edge, edge->opposite(v));
                    }
                }
            }
        };

        ///////////////////////
        // Graph Contraction //
        ///////////////////////

        template <>
        struct GraphContraction<Undirected>
        {
            template <typename Vertex, typename Edge>
            static std::pair<Vertex*, Edge*> splitEdge(Edge * e)
            {
                auto * v = e->graph()->addVertex(); 
                auto * k = e->graph()->addEdge(v, e->vertex2());
                e->setVertex2(v);
                return std::make_pair(v, k);
            }

            template <typename Vertex, typename Edge>
            static void contractEdge(Edge * e, Vertex * v1, Vertex * v2)
            {
                for (auto * k : boost::adaptors::reverse(v2->edges())) {
                    if (k != e) {
                        if (k->vertex1() == v2) {
                            k->setVertex1(v1);
                        } 
                        if (k->vertex2() == v2) {
                            k->setVertex2(v1);
                        }
                    }
                }
                e->graph()->removeVertex(v2);
            }

            template <typename Vertex>
            static Vertex * splitVertex(Vertex* v)
            {
                typedef GraphElementData<typename Vertex::DataType, typename Vertex::CustomDataType> VertexElementData;
                typedef GraphElementData<typename Vertex::Edge::DataType, typename Vertex::Edge::CustomDataType> EdgeElementData;
                Vertex * k = v->graph()->addVertex();
                static_cast<VertexElementData&>(*k) = static_cast<const VertexElementData&>(*v);
                for (auto * e : v->edges()) {
                    if (e->vertex1() == v && e->vertex2() != v) {
                        auto * f = v->graph()->addEdge(k, e->vertex2());
                        static_cast<EdgeElementData&>(*f) = static_cast<const EdgeElementData&>(*e);
                    } else if (e->vertex1() != v && e->vertex2() == v) {
                        auto * f = v->graph()->addEdge(e->vertex1(), k);
                        static_cast<EdgeElementData&>(*f) = static_cast<const EdgeElementData&>(*e);
                    } else {
                        auto * f = v->graph()->addEdge(k, k);
                        static_cast<EdgeElementData&>(*f) = static_cast<const EdgeElementData&>(*e);
                    }
                }
                return k;
            }
        };

        ////////////////////////////////////
        // Standard Edge Isomorphism Test //
        ////////////////////////////////////

        template <>
        struct StandardEdgeIsomorphismTest<Undirected>
        {
            template <typename Vertex, typename Edge>
            bool operator() (const Vertex * v1, const Vertex * v2, const Vertex * w1, const Vertex * w2, const Edge * e, const Edge * f) const
            {
                return ((v1 == e->vertex1() && v2 == e->vertex2()) || (v1 == e->vertex2() && v2 == e->vertex1())) 
                    && ((w1 == f->vertex1() && w2 == f->vertex2()) || (w1 == f->vertex2() && w2 == f->vertex1()));
            }
        };
    }
}

#endif // GULO_GRAPH_EDGEUNDIRECTED_HPP