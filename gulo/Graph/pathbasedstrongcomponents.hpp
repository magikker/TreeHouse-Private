#ifndef GULO_GRAPH_PATHBASEDSTRONGCOMPONENTS_HPP
#define GULO_GRAPH_PATHBASEDSTRONGCOMPONENTS_HPP

#include <stack>
#include <vector>
#include <set>

namespace Gulo
{
    namespace Detail
    {
        template <typename VertexPointer, typename EdgePointer>
        struct PathBasedStrongComponentGatherer : public DepthFirstVisitor
        {
            PathBasedStrongComponentGatherer(size_t i)
            {
                preorderNumber.resize(i, -1);
            }

            void reachUnvisitedVertex(VertexPointer v)
            {
                preorderNumber[v->id()] = c;
                ++c;
                P.push(v);
                S.push(v);
            }

            void finishVertex(VertexPointer v)
            {
                if (P.size() > 0 && P.top() == v) {
                    components.push_back({});
                    while (S.top() != v) {
                        components.back().insert(S.top());
                        S.pop();
                    }
                    components.back().insert(S.top());
                    S.pop();
                    P.pop();
                }
            }

            void reachQueuedVertex(VertexPointer v)
            {
                while (P.size() > 0 && preorderNumber[P.top()->id()] > preorderNumber[v->id()]) {
                    P.pop();
                }
            }

            size_t c{0};
            std::vector<int> preorderNumber;
            std::stack<VertexPointer> P;
            std::stack<VertexPointer> S;
            std::vector<std::set<VertexPointer> > components;
        };

        template <typename VertexPointer, typename EdgePointer>
        struct PathBasedStrongComponentTester : public DepthFirstVisitor
        {
            PathBasedStrongComponentTester(size_t i)
            {
                preorderNumber.resize(i, -1);
            }

            void reachUnvisitedVertex(VertexPointer v)
            {
                preorderNumber[v->id()] = c;
                ++c;
                P.push(v);
                S.push(v);
            }

            void finishVertex(VertexPointer v)
            {
                if (P.size() > 0 && P.top() == v) {
                    ++components;
                    while (S.top() != v) {
                        S.pop();
                    }
                    S.pop();
                    P.pop();
                }
            }

            void reachQueuedVertex(VertexPointer v)
            {
                while (P.size() > 0 && preorderNumber[P.top()->id()] > preorderNumber[v->id()]) {
                    P.pop();
                }
            }

            bool terminate() const 
            {
                return components > 1;
            }

            size_t c{0};
            std::vector<int> preorderNumber;
            std::stack<VertexPointer> P;
            std::stack<VertexPointer> S;
            size_t components{0};
        };
    }
}

#endif // GULO_GRAPH_PATHBASEDSTRONGCOMPONENTS_HPP
