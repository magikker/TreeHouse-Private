#ifndef GULO_PRIORITYQUEUE_HPP
#define GULO_PRIORITYQUEUE_HPP

#include <boost/integer/static_log2.hpp>
#include <type_traits>

// implements a priority queue as a 4-ary heap
// priorities are templated by type T
// the back of the queue is used to store a postorder list
// size_t is used as the index type

namespace Gulo
{
    template <bool IsConst, typename T, typename Comparator = std::less<T> >
    class PriorityQueue
    {            
        static const size_t D = 4;
        static const size_t Log2D = boost::static_log2<D>::value;

        
    public:
        
        PriorityQueue(T defaultValue, size_t elementCount)
        : values((T*)malloc(elementCount * sizeof(T))),
          data((size_t*)malloc((elementCount+1) * sizeof(size_t))),
          index((size_t*)malloc(elementCount * sizeof(size_t))),
          end(0),
          storage(elementCount+1)
        {
            for (size_t i = 0; i < elementCount; ++i) {
                values[i] = defaultValue;
            }
        }
        
        ~PriorityQueue()
        {
            free(values);
            free(data);
            free(index);
        }
              
        size_t size() const
        {
            return end;
        }
        
        bool empty() const
        {
            return end == 0;
        }
        
        void push(size_t id, T value)
        {
            values[id] = value;
            size_t pos = end;
            size_t up = parent(pos);
            ++end;
            while (pos != 0 && up < end && comp(value, values[data[up]])) {
                data[pos] = data[up];
                index[data[pos]] = pos;
                pos = up;
                up = parent(pos); 
            }
            data[pos] = id;
            index[id] = pos;
        }
        
        void increase(size_t id, T value)
        {
            if (values[id] != value) {
                values[id] = value;
                size_t pos = index[id];
                size_t up = parent(pos);
                while (pos != 0 && up < end && comp(value, values[data[up]])) {
                    data[pos] = data[up];
                    index[data[pos]] = pos;
                    pos = up;
                    up = parent(pos); 
                }
                data[pos] = id;
                index[id] = pos;
            }
        }
        
        size_t front() const
        {
            return data[0];
        }
        
        void pop()
        {
            --storage;
            data[storage] = data[0];
            --end;      
            size_t pos = end;
            size_t swap = 0;
            while (pos != swap) {
                std::swap(data[pos], data[swap]);
                std::swap(index[data[pos]], index[data[swap]]);
                pos = swap;
                for (int i = 0; i < D && child(pos)+i < end; ++i)  {
                    if (comp(values[data[child(pos) + i]],values[data[swap]])) {
                        swap = child(pos) + i;
                    }
                }
            };
        }
        
        T priority(size_t id) const
        {
            return values[id];
        }
        
        size_t at(size_t i) const
        {
            return data[i];
        }
        
        size_t storagePos() const
        {
            return storage;
        }
               
    private:
        
        size_t parent(size_t i) const
        {
            return ((i-1)>>Log2D);
        }
        
        size_t child(size_t i) const
        {
            return (i<<Log2D)+1;
        }
         
        T * values;          
        size_t * data;
        size_t * index;
        Comparator comp;
        size_t end;
        size_t storage;
        size_t maxElements;
    };
}

#endif // GULO_PRIORITYQUEUE_HPP