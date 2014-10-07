#ifndef GULO_DIJKSTRAVISITOR_HPP
#define GULO_DIJKSTRAVISITOR_HPP

#include <cassert>
#include "priorityqueue.hpp"
#include "breadthfirsttraversal.hpp"

namespace Gulo
{     
    struct Dijkstra 
    {
        template <bool IsConst, typename Visitor, typename Vertex, typename Edge>
        struct InternalVisitor
        {
            typedef typename std::conditional<IsConst, const Vertex *, Vertex *>::type VertexPointer;
            typedef typename std::conditional<IsConst, const Edge *, Edge *>::type EdgePointer;
            typedef typename Visitor::LengthType LengthType;
            
            InternalVisitor(Visitor & vis, size_t vertexCount)
            : visitor(vis),
              queue(0, vertexCount),
              u(0),
              len(0)
            {
               if (visitor.pathLengths.size() == 0) {
                   visitor.pathLengths = std::vector<LengthType>(vertexCount);
               }
               if (visitor.predecessors.size() == 0) {
                   visitor.predecessors = std::vector<const void*>(vertexCount, 0);
               }
            }
            
            ~InternalVisitor()
            {
            }
            
            void startTraversal(VertexPointer vertex)
            {
                u = vertex;
                len = 0;
                visitor.startTraversal(vertex);
            }
            
            void finishTraversal(VertexPointer vertex)
            {
                visitor.finishTraversal(vertex);
            }
            
            void startVertex(VertexPointer vertex)
            {
                u = vertex;
                visitor.startVertex(vertex);
            }
            
            void finishVertex(VertexPointer vertex)
            {
                visitor.finishVertex(vertex);
            }
            
            void startEdge(EdgePointer edge)
            {
                edgePrev = edge;
                visitor.startEdge(edge);
                len = visitor.getLength(edge);
                assert(len >= 0); // todo, throw exception if len < 0
            }
            
            void finishEdge(EdgePointer edge)
            {
                visitor.finishEdge(edge);
            }
            
            void reachUnvisitedVertex(VertexPointer vertex)
            {
                visitor.edgePred[vertex->id()] = (const void*)edgePrev;
                queue.push(vertex->id(), queue.priority(u->id()) + len);
                visitor.reachUnvisitedVertex(vertex);
            }
            
            void reachQueuedVertex(VertexPointer vertex)
            {
                if (queue.priority(u->id()) + len < queue.priority(vertex->id())) {
                    queue.increase(vertex->id(), queue.priority(u->id()) + len);
                    visitor.predecessors[vertex->id()] = u; 
                }
                visitor.reachQueuedVertex(vertex);
            }
            
            void reachVisitedVertex(VertexPointer vertex)
            {
                visitor.reachVisitedVertex(vertex);
            }
            
            void reachBlockedVertex(VertexPointer vertex)
            {
                visitor.reachBlockedVertex(vertex);
            }
            
            std::pair<VertexPointer, EdgePointer> queueFront() const
            {
                auto it = queue.front();
                visitor.pathLengths[it] = queue.priority(it);
                return std::make_pair(u->graph()->vertex(it), const_cast<EdgePointer>(reinterpret_cast<const Edge*>(visitor.edgePred[u->graph()->vertex(it)->id()])));
            }
            
            void queuePop()
            {
                queue.pop();
            }
            
            bool queueEmpty() const 
            {
                return queue.empty();
            }
            
            bool terminate() const
            {
                return visitor.terminate();
            }
            
            Visitor & visitor;
            PriorityQueue<IsConst, LengthType> queue;
            VertexPointer u;
            LengthType len;
            EdgePointer edgePrev;
        };
        
        template <typename DirectedOption, typename Visitor, typename Vertex>
        static void visit(Vertex * vertex, Visitor & visitor, char * visitorState = 0)
        {  
            InternalVisitor<false, Visitor, Vertex, typename Vertex::Graph::Edge> vis(visitor, vertex->graph()->vertices().size());
            if (visitor.edgePred.size() == 0) {
                visitor.edgePred.assign(vertex->graph()->vertices().size(), nullptr);
            }
            Detail::BreadthFirstTraversal::template go<DirectedOption>(vis, vertex, visitorState);
        }
        
        template <typename DirectedOption, typename Visitor, typename Vertex>
        static void visit(const Vertex * vertex, Visitor & visitor, char * visitorState = 0)
        {  
            InternalVisitor<true, Visitor, const Vertex, const typename Vertex::Graph::Edge> vis(visitor, vertex->graph()->vertices().size());
            if (visitor.edgePred.size() == 0) {
                visitor.edgePred.assign(vertex->graph()->vertices().size(), nullptr);
            }
            Detail::BreadthFirstTraversal::template go<DirectedOption>(vis, vertex, visitorState);
        }
    };
    
    template <typename LengthT = double>
    struct DijkstraVisitor
    {
        typedef LengthT LengthType;
        typedef Dijkstra Traversal;
        
        friend class Dijkstra;

        template <typename Vertex> void startTraversal(Vertex * vertex) {}
        template <typename Vertex> void finishTraversal(Vertex * vertex) {}
        template <typename Vertex> void startVertex(Vertex * vertex) {}
        template <typename Vertex> void finishVertex(Vertex * vertex) {}
        template <typename Edge>   void startEdge(Edge * edge) {}
        template <typename Edge>   void finishEdge(Edge * edge) {}
        template <typename Vertex> void reachUnvisitedVertex(Vertex * vertex) {}
        template <typename Vertex> void reachVisitedVertex(Vertex * vertex) {}
        template <typename Vertex> void reachQueuedVertex(Vertex * vertex) {}
        template <typename Vertex> void reachBlockedVertex(Vertex * vertex) {}
        
        template <typename Vertex> LengthType pathLength(const Vertex * vertex) {return pathLengths[vertex->id()];}
        template <typename Vertex> Vertex * predecessor(Vertex * v) const {return const_cast<Vertex*>(reinterpret_cast<const Vertex*>(predecessors[v->id()]));}
        template <typename Vertex> const Vertex * predecessor(const Vertex * v) const {return reinterpret_cast<const Vertex*>(predecessors[v->id()]);}
        
        bool terminate() const {return false;}

    protected:
        
        std::vector<LengthType> pathLengths;
        std::vector<const void*> predecessors;  // dijkstra predecessors
        std::vector<const void*> edgePred;      // predecessor edges in the breadth first traversal (prior to queueing)
    };
    
    /* example visitors
        
        typedef MyGraph::Edge Edge;
        typedef MyGraph::Vertex Vertex;

        struct Visitor : public DijkstraVisitor<double>
        {
            double getLength(const Edge * e) {return e->length();}; <- must provide a getLength function for your graph

            void reachUnvisitedVertex(const Vertex * vertex) {std::cout << vertex->value() << std::endl;}
        };


        ----------------


        template<typename Length, typename Graph>
        struct MoreGenericVisitor : public DijkstraVisitor<Length>
        {
            typedef typename Graph::Vertex Vertex;
            typedef typename Graph::Edge Edge;

            MoreGenericVisitor(const std::function<Length(const Edge *)> & functor)
            : getLength(functor)
            {
            }

            void reachUnvisitedVertex(const Vertex * vertex) {std::cout << vertex->value() << std::endl;}

            const std::function<Length(const Edge *)> & getLength;
        };

        usage of more generic version:
        MoreGenericVisitor<double, MyGraph> visitor([](const Edge* edge)->double{return edge->length();};);

    */
}

#endif // GULO_DIJKSTRAVISITOR_HPP
