#ifndef GULO_GRAPH_ELEMENTVECTORACCESSOR_HPP
#define GULO_GRAPH_ELEMENTVECTORACCESSOR_HPP

#include "elementvector.hpp"

namespace Gulo
{
    namespace Detail
    {
        struct GraphElementVectorAccessor
        {
            template <typename T, typename Alloc>
            static std::vector<T*, Alloc> & getVector(GraphElementVector<T, Alloc> & vector)
            {
                return vector.vector();
            }
        };
    }
}

#endif // GULO_GRAPH_ELEMENTVECTORACCESSOR_HPP