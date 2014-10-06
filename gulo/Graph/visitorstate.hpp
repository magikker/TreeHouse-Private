#ifndef GULO_GRAPH_VISITORSTATE_HPP
#define GULO_GRAPH_VISITORSTATE_HPP

namespace Gulo
{
    namespace Detail
    {
        // State of a vertex or edge during traversal of a visitor over a graph 
        struct VisitorState
        {
            static const char Unvisited = 0;
            static const char Queued = 1;
            static const char Visited = 2;
            static const char Blocked = 3;
        };
    }
}

#endif // GULO_GRAPH_VISITORSTATE_HPP
