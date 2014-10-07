#ifndef GULO_GRAPH_OUTDIRECTED_HPP
#define GULO_GRAPH_OUTDIRECTED_HPP

#include "elementvectoraccessor.hpp"
#include "elementvector.hpp"
#include "forward.hpp"
#include "graphaccessor.hpp"
#include "elementaccessor.hpp"
#include "graphtraversalout.hpp"
#include "elementdata.hpp"
#include "weaklyconnectedcomponentsoutdirected.hpp"
#include "pathbasedstrongcomponents.hpp"
#include "arborescencetest.hpp"

#include <algorithm>

namespace Gulo
{
    template <typename StoragePolicy> struct Directed;

    namespace Detail
    { 
        /////////////////////
        // Graph Extension //
        /////////////////////

        template <typename Graph, typename Vertex>
        class GraphExtension<Graph, Vertex, Directed<StoreOutEdges>>
        {
        public:
            
            constexpr static bool isUndirectedGraph() {return false;}
            constexpr static bool isDirectedGraph()   {return true;}
            constexpr static bool isBidirectedGraph() {return false;}
            constexpr static bool storesOutEdges()    {return true;}
            constexpr static bool storesInEdges()     {return false;} 

            std::vector<std::set<Vertex*>> weaklyConnectedComponents()
            {
                OutDirectedWeaklyConnectedComponentGatherer<Vertex*> v;
                static_cast<Graph*>(this)->visit(v);
                v.components.erase(std::remove_if(v.components.begin(), 
                              v.components.end(),
                              [](const std::set<Vertex*> & x){return x.empty();}),
                           v.components.end());
                return v.components;
            }

            std::vector<std::set<const Vertex*>> weaklyConnectedComponents() const
            {
                OutDirectedWeaklyConnectedComponentGatherer<const Vertex*> v;
                static_cast<const Graph*>(this)->visit(v);
                v.components.erase(std::remove_if(v.components.begin(), 
                              v.components.end(),
                              [](const std::set<const Vertex*> & x){return x.empty();}),
                           v.components.end());
                return v.components;
            }

            bool isWeaklyConnectedGraph() const
            {
                return weaklyConnectedComponents().size() == 1;
            }

            std::vector<std::set<Vertex*>> stronglyConnectedComponents()
            {
                PathBasedStrongComponentGatherer<Vertex*, typename Vertex::Edge*> v(static_cast<Graph*>(this)->vertices().size());
                static_cast<Graph*>(this)->visit(v);
                return v.components;
            }

            std::vector<std::set<const Vertex*>> stronglyConnectedComponents() const
            {
                PathBasedStrongComponentGatherer<const Vertex*, const typename Vertex::Edge*> v(static_cast<const Graph*>(this)->vertices().size());
                static_cast<const Graph*>(this)->visit(v);
                return v.components;
            } 

            bool isStronglyConnectedGraph() const
            {
                PathBasedStrongComponentTester<const Vertex*, const typename Vertex::Edge*> v(static_cast<const Graph*>(this)->vertices().size());
                static_cast<const Graph*>(this)->visit(v);
                return !(v.terminate());
            }

            bool isArborescence() const
            {
                return (arborescenceRoot() != 0);
            }

            Vertex * arborescenceRoot()
            {
                return const_cast<Vertex*>(const_cast<const Graph*>(this)->arborescenceRoot());
            }

            const Vertex * arborescenceRoot() const
            {   
                size_t len = static_cast<const Graph*>(this)->vertices().size();
                char * visitorState = (char*)calloc(len, sizeof(char)); 
                for (auto * e : static_cast<const Graph*>(this)->edges()) {
                    if (e != nullptr) {
                        visitorState[e->vertex2()->id()] = 1;
                    }
                }
                const Vertex * root = nullptr;
                for (size_t i = 0; i < len; ++i) {
                    if (visitorState[i] == 0) {
                        root = static_cast<const Graph*>(this)->vertices()[i];
                        break;
                    }
                }
                if (root == nullptr) {
                    free(visitorState);
                    return nullptr;
                } 
                memset (visitorState, 0, len);
                ArborescenceTester<const Vertex*> vis;
                ArborescenceTester<const Vertex*>::Traversal::template visit<Directed<StoreOutEdges>>(root, vis, visitorState);
                for (size_t i = 0; i < len; ++i) {
                    if (visitorState[i] == 0 && static_cast<const Graph*>(this)->vertices()[i] != 0) {
                        free(visitorState);
                        return nullptr;
                    }
                }
                free(visitorState);
                return (vis.terminate()) ? nullptr : root;
            }  

            bool isTree() const
            {
                return static_cast<const Graph*>(this)->isAcyclicGraph() && static_cast<const Graph*>(this)->isWeaklyConnectedGraph();
            }
        };

        ////////////////////
        // Edge Extension //
        ////////////////////

        template <typename Vertex, typename Edge>
        class EdgeExtension<Vertex, Edge, Directed<StoreOutEdges> > : protected GraphElementVectorAccessor
        {
        public:

            void flipDirection()
            {
                if (!static_cast<Edge*>(this)->isLoop()) {
                    auto & outVec = getVector(static_cast<Edge*>(this)->m_vertex1->outgoingEdges());
                    outVec.erase(std::find(outVec.begin(), outVec.end(), this));
                    getVector(static_cast<Edge*>(this)->m_vertex2->outgoingEdges()).push_back(static_cast<Edge*>(this));
                    std::swap(static_cast<Edge*>(this)->m_vertex1, static_cast<Edge*>(this)->m_vertex2);
                }
            }
        };

        /////////////////////////
        // Copy Edge Extension //
        /////////////////////////

        template <>
        struct CopyEdgeExtension<Directed<StoreOutEdges> >
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
        struct CopyVertexExtension<Directed<StoreOutEdges>>
        {
            template <typename Edge>
            static void copy(Edge * target, const Edge * source)
            {
            }
        };

        //////////////////////
        // Vertex Extension //
        //////////////////////

        template <typename Edge>
        class VertexExtension<Edge, Directed<StoreOutEdges> > 
        {
            typedef typename Edge::Vertex Vertex;

        public:

            GraphElementVector<Edge> & outgoingEdges()             {return m_edges;}
            const GraphElementVector<Edge> & outgoingEdges() const {return m_edges;}

            std::vector<const Edge *> edgesWith(const Vertex * v) const
            {
                std::vector<const Edge *> result;
                if (v != 0) {
                    std::for_each(m_edges.begin(), m_edges.end(), [&](const Edge * e)->void{ if (!e->isLoop() && v == e->vertex2()) result.push_back(e);});
                    std::for_each(v->m_edges.begin(), v->m_edges.end(), [&](const Edge * e)->void{ if (this == e->vertex2()) result.push_back(e);});                  
                }
                return result;
            }
            
            std::vector<Edge *> edgesWith(const Vertex * v)
            {
                std::vector<Edge *> result;
                if (v != 0) {
                    std::for_each(m_edges.begin(), m_edges.end(), [&](const Edge * e)->void{ if (!e->isLoop() && v == e->vertex2()) result.push_back(const_cast<Edge*>(e));});
                    std::for_each(v->m_edges.begin(), v->m_edges.end(), [&](const Edge * e)->void{ if (this == e->vertex2()) result.push_back(const_cast<Edge*>(e));});                  
                }
                return result;
            }

            bool isOutAdjacentTo(const Vertex * other) const
            {
                for (auto * e : m_edges) {
                    if (other == e->vertex2()) {
                       return true;
                    }
                }
                return false;
            }

            bool isInAdjacentTo(const Vertex * other) const
            {
                for (auto * e : other->m_edges) {
                    if (this == e->vertex2()) {
                       return true;
                    }
                }
                return false;
            }

            bool isAdjacentTo(const Vertex * other) const
            {
                return isInAdjacentTo(other) || isOutAdjacentTo(other);
            }

            size_t outDegree() const
            {
                return m_edges.size();
            }

        private:

            GraphElementVector<Edge> m_edges;
        };

        ///////////////////////////
        // Modify Vertex Storage //
        ///////////////////////////

        template <>
        struct ModifyVertexStorage<Directed<StoreOutEdges>> : protected GraphElementVectorAccessor,
                                                              protected GraphAccessor,
                                                              protected GraphElementAccessor
        {
            template <typename Edge>
            static void addEdgeToVertexStorage(Edge * e) 
            {
                getVector(e->vertex1()->outgoingEdges()).push_back(e);
            }

            template <typename Edge>
            static void removeEdgeFromVertexStorage(Edge * e) 
            {
                if (e->vertex1() != nullptr) {
                    auto & vec = getVector(e->vertex1()->outgoingEdges());
                    vec.erase(std::find(vec.begin(), vec.end(), e));
                }
            }

            template <typename Edge, typename Vertex>
            static void setVertex1Storage(Edge * e, Vertex * v)
            {
                auto & vec = getVector(e->vertex1()->outgoingEdges());
                vec.erase(std::find(vec.begin(), vec.end(), e));
                getVector(v->outgoingEdges()).push_back(e);
            }

            template <typename Edge, typename Vertex>
            static void setVertex2Storage(Edge * e, Vertex * v)
            {
            }

            template <typename Vertex>
            static void disconnect(Vertex * v)
            {
                typedef typename Vertex::Edge Edge;
                auto & edgeVec = getVector(v->graph()->edges());
                bool modifiedEdgeList = false;
                for (size_t i = 0; i < edgeVec.size(); ++i) {
                    if (edgeVec[i] != 0 && (edgeVec[i]->vertex1() == v || edgeVec[i]->vertex2() == v)) {
                        if (edgeVec[i]->vertex2() == v) {
                            auto & E = getVector(edgeVec[i]->vertex1()->outgoingEdges());
                            E.erase(std::find(E.begin(), E.end(), edgeVec[i]));
                        }
                        edgeVec[i]->~Edge();
                        edgePool(v->graph())->free(edgeVec[i]);
                        if (!v->graph()->elementReindexingIsSuspended()) {
                            edgeVec.erase(edgeVec.begin() + i);
                            --i;
                            modifiedEdgeList = true;
                        } else {
                            *(edgeVec.begin() + i) = 0;
                        }
                    } else if (modifiedEdgeList && !v->graph()->elementReindexingIsSuspended()) {
                        getElementId(edgeVec[i]) = i;
                    }
                }
                getVector(v->outgoingEdges()).clear();
            }
        };

        /////////////////////////////
        // Parallel Edge Detection //
        /////////////////////////////

        template <typename Edge>
        struct ParallelEdgeDetection<Edge, Directed<StoreOutEdges> >
        {
            bool existsCoParallelEdge() const
            {
                for (auto * e : static_cast<const Edge*>(this)->vertex1()->outgoingEdges()) {
                    if (e != this && (e->vertex1() == static_cast<const Edge*>(this)->vertex1() && e->vertex2() == static_cast<const Edge*>(this)->vertex2())) {
                        return true;
                    }
                }
                return false;
            }

            bool existsCounterParallelEdge() const
            {
                for (auto * e : static_cast<const Edge*>(this)->vertex2()->outgoingEdges()) {
                    if (e->vertex1() == static_cast<const Edge*>(this)->vertex2() && e->vertex2() == static_cast<const Edge*>(this)->vertex1()) {
                        return true;
                    }
                }
                return false;
            }

            bool existsParallelEdge() const
            {
                return (existsCoParallelEdge() || existsCounterParallelEdge());
            }
        };

        ///////////////////////
        // Graph Contraction //
        ///////////////////////

        template <>
        struct GraphContraction<Directed<StoreOutEdges>>
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
               for (auto * k : e->graph()->edges()) {
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
                Vertex * q = v->graph()->addVertex();
                static_cast<VertexElementData&>(*q) = static_cast<const VertexElementData&>(*v);
                size_t k = v->graph()->edges().size();
                for (size_t i = 0; i < k; ++i) {
                    auto * e = v->graph()->edges()[i];
                    if (e != nullptr) {
                        if (e->vertex1() == v && e->vertex2() == v) {
                            auto * f = v->graph()->addEdge(q, q);
                            static_cast<EdgeElementData&>(*f) = static_cast<const EdgeElementData&>(*e);
                        } else if (e->vertex1() == v) {
                            auto * f = v->graph()->addEdge(q, e->vertex2());
                            static_cast<EdgeElementData&>(*f) = static_cast<const EdgeElementData&>(*e);
                        } else if (e->vertex2() == v) {
                            auto * f = v->graph()->addEdge(e->vertex1(), q);
                            static_cast<EdgeElementData&>(*f) = static_cast<const EdgeElementData&>(*e);
                        }
                    }
                }
                return q;
            }
        };

        ////////////////////////////////////
        // Standard Edge Isomorphism Test //
        ////////////////////////////////////

        template <>
        struct StandardEdgeIsomorphismTest<Directed<StoreOutEdges>>
        {
            template <typename Vertex, typename Edge>
            bool operator() (const Vertex * v1, const Vertex * v2, const Vertex * w1, const Vertex * w2, const Edge * e, const Edge * f) const
            {
                return (v1 == e->vertex1() && v2 == e->vertex2() && w1 == f->vertex1() && w2 == f->vertex2()) 
                    || (v1 == e->vertex2() && v2 == e->vertex1() && w1 == f->vertex2() && w2 == f->vertex1());
            }
        };
    }
}

#endif // GULO_GRAPH_OUTDIRECTED_HPP