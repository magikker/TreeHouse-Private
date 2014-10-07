#ifndef GULO_GRAPH_INOUTBIDIRECTED_HPP
#define GULO_GRAPH_INOUTBIDIRECTED_HPP

#include "elementvectoraccessor.hpp"
#include "elementvector.hpp"
#include "forward.hpp"
#include "graphaccessor.hpp"
#include "elementaccessor.hpp"
#include "exception.hpp"
#include "graphtraversal.hpp"
#include "graphtraversalinout.hpp"
#include "elementdata.hpp"
#include <boost/range/adaptor/reversed.hpp>
#include "weaklyconnectedcomponentsundirected.hpp"
#include "pathbasedstrongcomponents.hpp"
#include "arborescencetest.hpp"

#include <algorithm>

namespace Gulo
{
    template <typename StoragePolicy> struct Bidirected;

    namespace Detail
    { 
        /////////////////////
        // Graph Extension //
        /////////////////////

        template <typename Graph, typename Vertex>
        class GraphExtension<Graph, Vertex, Bidirected<StoreInOutEdges>>
        {
        public:

            template <typename Visitor>
            void visitUndirected(Visitor && visitor)
            {
                GraphTraversal<InOutUndirected<Bidirected<StoreInOutEdges>>>::template visit<Graph>(static_cast<Graph&>(*this), visitor);
            }

            template <typename Visitor>
            void visitUndirected(Visitor && visitor) const
            {
                GraphTraversal<InOutUndirected<Bidirected<StoreInOutEdges>>>::template visit<Graph>(static_cast<const Graph&>(*this), visitor);
            }

            constexpr static bool isUndirectedGraph() {return false;}
            constexpr static bool isDirectedGraph()   {return false;}
            constexpr static bool isBidirectedGraph() {return true;}
            constexpr static bool storesOutEdges()    {return true;}
            constexpr static bool storesInEdges()     {return true;} 

            std::vector<std::set<Vertex*>> weaklyConnectedComponents()
            {
                UndirectedWeaklyConnectedComponentGatherer<Vertex*> v;
                visitUndirected(v);
                return v.components;
            }

            std::vector<std::set<const Vertex*>> weaklyConnectedComponents() const
            {
                UndirectedWeaklyConnectedComponentGatherer<const Vertex*> v;
                visitUndirected(v);
                return v.components;
            } 

            bool isWeaklyConnectedGraph() const
            {
                UndirectedWeaklyConnectedGraphTester<const Vertex*> v;
                static_cast<const Graph*>(this)->visitUndirected(v);
                return !(v.terminate());
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
                ArborescenceTester<const Vertex*> vis;
                visit(vis);
                return (vis.terminate()) ? nullptr : vis.initialVertex;
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
        class EdgeExtension<Vertex, Edge, Bidirected<StoreInOutEdges> > : protected GraphElementVectorAccessor
        {
            friend struct ModifyVertexStorage<Bidirected<StoreInOutEdges>>;

        public:

            void flipDirection()
            {
                if (!isBidirected()) {
                    auto & outVec = getVector(static_cast<Edge*>(this)->vertex1()->outgoingEdges());
                    auto & inVec = getVector(static_cast<Edge*>(this)->vertex2()->incomingEdges());
                    inVec.erase(std::find(inVec.begin(), inVec.end(), this));
                    outVec.erase(std::find(outVec.begin(), outVec.end(), this));
                    getVector(static_cast<Edge*>(this)->vertex1()->incomingEdges()).push_back(static_cast<Edge*>(this));
                    getVector(static_cast<Edge*>(this)->vertex2()->outgoingEdges()).push_back(static_cast<Edge*>(this));
                    std::swap(static_cast<Edge*>(this)->m_vertex1, static_cast<Edge*>(this)->m_vertex2);
                }
            }

            bool isBidirected() const   {return m_isBidirected || static_cast<const Edge*>(this)->isLoop();}

            void setBidirected()
            {
                if (!isBidirected()) {
                    getVector(static_cast<Edge*>(this)->vertex1()->incomingEdges()).push_back(static_cast<Edge*>(this));
                    getVector(static_cast<Edge*>(this)->vertex2()->outgoingEdges()).push_back(static_cast<Edge*>(this));
                    m_isBidirected = true;
                }
            }

            void setUnidirected(Vertex * target)
            {
                if (isBidirected()) {
                    if (static_cast<Edge*>(this)->vertex1() == target) {
                        auto & outVec = getVector(static_cast<Edge*>(this)->vertex1()->outgoingEdges());
                        auto & inVec = getVector(static_cast<Edge*>(this)->vertex2()->incomingEdges());
                        inVec.erase(std::find(inVec.begin(), inVec.end(), this));
                        outVec.erase(std::find(outVec.begin(), outVec.end(), this));
                        std::swap(static_cast<Edge*>(this)->m_vertex1, static_cast<Edge*>(this)->m_vertex2);
                    } else if (static_cast<Edge*>(this)->vertex2() == target) {
                        auto & outVec = getVector(static_cast<Edge*>(this)->vertex2()->outgoingEdges());
                        auto & inVec = getVector(static_cast<Edge*>(this)->vertex1()->incomingEdges());
                        inVec.erase(std::find(inVec.begin(), inVec.end(), this));
                        outVec.erase(std::find(outVec.begin(), outVec.end(), this));
                    } else {
                        throw BadGraphElementException("Attempt to target a bidirected edge toward a vertex that is not a member of the edge");
                    }
                    m_isBidirected = false;
                } else if (!static_cast<Edge*>(this)->isBidirected() && target != static_cast<Edge*>(this)->vertex2()) {
                    flipDirection();
                }
            }

        private:

            bool m_isBidirected{false};
        };

        /////////////////////////
        // Copy Edge Extension //
        /////////////////////////

        template <>
        struct CopyEdgeExtension<Bidirected<StoreInOutEdges>>
        {
            template <typename Edge>
            static void copy(Edge * target, const Edge * source)
            {
                if (!source->isBidirected()) {
                    target->setUnidirected(target->vertex2());
                } else {
                    target->setBidirected();
                }
            }
        };

        ///////////////////////////
        // Copy Vertex Extension //
        ///////////////////////////

        template <>
        struct CopyVertexExtension<Bidirected<StoreInOutEdges>>
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
        class VertexExtension<Edge, Bidirected<StoreInOutEdges> > 
        {
            typedef typename Edge::Vertex Vertex;

        public:

            template <typename Visitor>
            void visitUndirected(Visitor && visitor)          {std::remove_reference<Visitor>::type::Traversal::template visit<InOutUndirected<Bidirected<StoreInOutEdges>>>(static_cast<typename Edge::Vertex*>(this), visitor);}

            template <typename Visitor>
            void visitUndirected(Visitor && visitor) const    {std::remove_reference<Visitor>::type::Traversal::template visit<InOutUndirected<Bidirected<StoreInOutEdges>>>(static_cast<const typename Edge::Vertex*>(this), visitor);}

            
            GraphElementVector<Edge> & incomingEdges()             {return m_in_edges;}
            const GraphElementVector<Edge> & incomingEdges() const {return m_in_edges;}

            GraphElementVector<Edge> & outgoingEdges()             {return m_out_edges;}
            const GraphElementVector<Edge> & outgoingEdges() const {return m_out_edges;}

            std::vector<const Edge *> edgesWith(const Vertex * v) const
            {
                std::vector<const Edge *> result;
                if (v != 0) {
                    std::for_each(this->m_out_edges.begin(), this->m_out_edges.end(), [&](const Edge*e)->void{ if (!e->isBidirected() && v == e->vertex2()) result.push_back(e);});
                    std::for_each(v->m_out_edges.begin(), v->m_out_edges.end(), [&](const Edge*e)->void{ if ( this == e->vertex2()) result.push_back(e);});                  
                }
                return result;
            }
            
            std::vector<Edge *> edgesWith(const Vertex * v)
            {
                std::vector<Edge *> result;
                if (v != 0) {
                    std::for_each(this->m_out_edges.begin(), this->m_out_edges.end(), [&](const Edge*e)->void{ if (!e->isBidirected() && v == e->vertex2()) result.push_back(const_cast<Edge*>(e));});
                    std::for_each(v->m_out_edges.begin(), v->m_out_edges.end(), [&](const Edge*e)->void{ if (this == e->vertex2()) result.push_back(const_cast<Edge*>(e));});                  
                }
                return result;
            }

            bool isOutAdjacentTo(const Vertex * other) const 
            {
                for (auto * e : this->m_out_edges) {
                    if (other == e->vertex2()) {
                       return true;
                    }
                }
                return false;
            }

            bool isInAdjacentTo(const Vertex * other) const
            {
                return (other->isOutAdjacentTo(static_cast<const Vertex*>(this)));
            }

            bool isBidirectionallyAdjacentTo(const Vertex * other) const
            {
                for (auto * e : other->m_out_edges) {
                    if ((this == e->vertex2() || this == e->vertex1()) && e->isBidirected()) {
                       return true;
                    }
                }
                return false;
            }

            bool isAdjacentTo(const Vertex * other) const
            {
                return isInAdjacentTo(other) || isOutAdjacentTo(other);
            }
            
            size_t inDegree() const
            {
                return m_in_edges.size();
            }

            size_t outDegree() const
            {
                return m_out_edges.size();
            }

            size_t degree() const
            {
                return m_in_edges.size() + m_out_edges.size();
            }

        private:

            GraphElementVector<Edge> m_in_edges;
            GraphElementVector<Edge> m_out_edges;
        };

        ///////////////////////////
        // Modify Vertex Storage //
        ///////////////////////////

        template <>
        struct ModifyVertexStorage<Bidirected<StoreInOutEdges>> : protected GraphElementVectorAccessor,
                                                                  protected GraphAccessor,
                                                                  protected GraphElementAccessor
        {
            template <typename Edge>
            static void addEdgeToVertexStorage(Edge * e) 
            {
                getVector(e->vertex1()->outgoingEdges()).push_back(e);
                getVector(e->vertex2()->incomingEdges()).push_back(e);
            }

            template <typename Edge>
            static void removeEdgeFromVertexStorage(Edge * e) 
            {
                if (!e->isLoop() && e->isBidirected()) {
                    if (e->vertex2() != nullptr) {
                        auto & vec = getVector(e->vertex2()->outgoingEdges());
                        vec.erase(std::find(vec.begin(), vec.end(), e));
                    }
                    if (e->vertex1() != nullptr) {
                        auto & vec = getVector(e->vertex1()->incomingEdges());
                        vec.erase(std::find(vec.begin(), vec.end(), e));
                    }
                }
                if (e->vertex1() != nullptr) {
                    auto & vec = getVector(e->vertex1()->outgoingEdges());
                    vec.erase(std::find(vec.begin(), vec.end(), e));
                }
                if (e->vertex2() != nullptr) {
                    auto & vec = getVector(e->vertex2()->incomingEdges());
                    vec.erase(std::find(vec.begin(), vec.end(), e));
                }
            }

            template <typename Edge, typename Vertex>
            static void setVertex1Storage(Edge * e, Vertex * v)
            {
                if (e->isBidirected() && !e->isLoop()) {
                    auto & vec = getVector(e->vertex1()->incomingEdges());
                    vec.erase(std::find(vec.begin(), vec.end(), e));
                    if (e->vertex2() != v) {
                        getVector(v->incomingEdges()).push_back(e);
                    }
                }
                auto & vec = getVector(e->vertex1()->outgoingEdges());
                vec.erase(std::find(vec.begin(), vec.end(), e));
                if (!(e->isBidirected() && e->vertex2() == v)) {
                    getVector(v->outgoingEdges()).push_back(e);
                }
                if (v == e->vertex2()) {
                    e->m_isBidirected = false;
                }
            }

            template <typename Edge, typename Vertex>
            static void setVertex2Storage(Edge * e, Vertex * v)
            {
                if (e->isBidirected() && !e->isLoop()) {
                    auto & vec = getVector(e->vertex2()->outgoingEdges());
                    vec.erase(std::find(vec.begin(), vec.end(), e));
                    if (e->vertex1() != v) {
                        getVector(v->outgoingEdges()).push_back(e);
                    }
                }
                auto & vec = getVector(e->vertex2()->incomingEdges());
                vec.erase(std::find(vec.begin(), vec.end(), e));
                if (!(e->isBidirected() && e->vertex1() == v)) {
                    getVector(v->incomingEdges()).push_back(e);
                }
                if (e->vertex1() == v) {
                    e->m_isBidirected = false;
                }
            }

            template <typename Vertex>
            static void disconnect(Vertex * v)
            {
                typedef typename Vertex::Edge Edge;
                auto & edgeVec = getVector(v->graph()->edges());
                size_t minEdgeIndex = v->graph()->edges().size();
                for (auto * e : v->incomingEdges()) {
                        minEdgeIndex = std::min(e->id(), minEdgeIndex);
                        edgeVec[e->id()] = 0;
                        if (e->vertex2() == v) {
                            auto & E = getVector(e->vertex1()->outgoingEdges());
                            E.erase(std::find(E.begin(), E.end(), e));
                            if (e->isBidirected()) {
                                auto & E2 = getVector(e->vertex1()->incomingEdges());
                                E2.erase(std::find(E2.begin(), E2.end(), e));
                            }
                            e->~Edge();
                            edgePool(v->graph())->free(e);
                        }
                }
                for (auto * e : v->outgoingEdges()) {
                    minEdgeIndex = std::min(e->id(), minEdgeIndex);
                    edgeVec[e->id()] = 0;
                    if (e->vertex1() == v) {
                            auto & E = getVector(e->vertex2()->incomingEdges());
                            E.erase(std::find(E.begin(), E.end(), e));
                            if (e->isBidirected()) {
                                auto & E2 = getVector(e->vertex2()->outgoingEdges());
                                E2.erase(std::find(E2.begin(), E2.end(), e));
                            }
                            e->~Edge();
                            edgePool(v->graph())->free(e);
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
                getVector(v->outgoingEdges()).clear();
                getVector(v->incomingEdges()).clear();
            }
        };

        /////////////////////////////
        // Parallel Edge Detection //
        /////////////////////////////

        template <typename Edge>
        struct ParallelEdgeDetection<Edge, Bidirected<StoreInOutEdges> >
        {
            bool existsCoParallelEdge() const
            {
                for (auto * e : static_cast<const Edge*>(this)->vertex1()->outgoingEdges()) {
                    if (e != this && ((e->vertex1() == static_cast<const Edge*>(this)->vertex1() && e->vertex2() == static_cast<const Edge*>(this)->vertex2()) || (e->isBidirected() && e->vertex1() == static_cast<const Edge*>(this)->vertex2() && e->vertex2() == static_cast<const Edge*>(this)->vertex1()))) {
                        return true;
                    }
                }
                return false;
            }

            bool existsCounterParallelEdge() const
            {
                for (auto * e : static_cast<const Edge*>(this)->vertex2()->outgoingEdges()) {
                    if (e != this && ((e->vertex1() == static_cast<const Edge*>(this)->vertex2() && e->vertex2() == static_cast<const Edge*>(this)->vertex1()) || (e->isBidirected() && e->vertex1() == static_cast<const Edge*>(this)->vertex1() && e->vertex2() == static_cast<const Edge*>(this)->vertex2()))) {
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
        struct GraphContraction<Bidirected<StoreInOutEdges>>
        {
            template <typename Vertex, typename Edge>
            static std::pair<Vertex*, Edge*> splitEdge(Edge * e)
            {
                auto * v = e->graph()->addVertex(); 
                auto * k = e->graph()->addEdge(v, e->vertex2());
                if (e->isBidirected()) {
                    k->setBidirected();
                }
                e->setVertex2(v);
                return std::make_pair(v, k);
            }

            template <typename Vertex, typename Edge>
            static void contractEdge(Edge * e, Vertex * v1, Vertex * v2)
            {
               for (auto * k : boost::adaptors::reverse(v2->outgoingEdges())) {
                    if (k != e) {
                        k->setVertex1(v1);
                    }
                }
                for (auto * k : boost::adaptors::reverse(v2->incomingEdges())) {
                    if (k != e) {
                        k->setVertex2(v1);
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
                for (auto * e : v->incomingEdges()) {
                    if (e->isLoop()) {
                        auto * f = v->graph()->addEdge(q, q);
                        static_cast<EdgeElementData&>(*f) = static_cast<const EdgeElementData&>(*e);
                    } else {
                        if (e->vertex2() == v) {
                            auto * f = v->graph()->addEdge(e->vertex1(), q);
                            static_cast<EdgeElementData&>(*f) = static_cast<const EdgeElementData&>(*e);
                            if (e->isBidirected()) f->setBidirected();
                        } else {
                            auto * f = v->graph()->addEdge(q, e->vertex2());
                            static_cast<EdgeElementData&>(*f) = static_cast<const EdgeElementData&>(*e);
                            f->setBidirected();
                        }
                    }
                }
                for (auto * e : v->outgoingEdges()) {
                    if (!e->isBidirected()) {
                        auto * f = v->graph()->addEdge(q, e->vertex2());
                        static_cast<EdgeElementData&>(*f) = static_cast<const EdgeElementData&>(*e);
                    }
                }
                return q;
            }
        };

        ////////////////////////////////////
        // Standard Edge Isomorphism Test //
        ////////////////////////////////////

        template <>
        struct StandardEdgeIsomorphismTest<Bidirected<StoreInOutEdges>>
        {
            template <typename Vertex, typename Edge>
            bool operator() (const Vertex * v1, const Vertex * v2, const Vertex * w1, const Vertex * w2, const Edge * e, const Edge * f) const
            {
                return (!e->isBidirected() && !f->isBidirected() && v1 == e->vertex1() && v2 == e->vertex2() && w1 == f->vertex1() && w2 == f->vertex2()) 
                      || (e->isBidirected() && f->isBidirected() && ((v1 == e->vertex1() && v2 == e->vertex2()) || (v1 == e->vertex2() && v2 == e->vertex1())) && ((w1 == f->vertex1() && w2 == f->vertex2()) || (w1 == f->vertex2() && w2 == f->vertex1())));
            }
        };
    }
}

#endif // GULO_GRAPH_INOUTDIRECTED_HPP