#ifndef GULO_GRAPH_ELEMENTVECTOR_HPP
#define GULO_GRAPH_ELEMENTVECTOR_HPP

#include <vector>

namespace Gulo
{
   template <typename T, typename Alloc = std::allocator<T*> >
   class GraphElementVector;
        
    namespace Detail
    {       
        template <bool IsConst, bool IsReverse, typename T, typename Alloc>
        class GraphElementVectorIterator
        {
            typedef typename std::conditional<IsReverse, typename std::vector<T*, Alloc>::const_reverse_iterator, typename std::vector<T*, Alloc>::const_iterator>::type Iterator;
            template <bool, bool, typename, typename> friend class GraphElementVectorIterator;
            friend class GraphElementVector<T, Alloc>;
            
        public: 
            
            typedef std::random_access_iterator_tag                             iterator_category;
            typedef typename std::conditional<IsConst, const T*, T*>::type      value_type;
            typedef typename std::vector<const T*>::difference_type             difference_type;
            typedef typename std::conditional<IsConst, const T*, T*>::type      pointer;
            typedef typename std::conditional<IsConst, const T*, T*>::type      reference;
            
            typedef typename std::remove_reference<typename std::vector<const T*, Alloc>::const_reference>::type reference_type;
            
            explicit GraphElementVectorIterator(Iterator it)
            : m_it(it)
            {
            }
            
            GraphElementVectorIterator(const GraphElementVectorIterator & other)
            : m_it(other.m_it)
            {
            }

            GraphElementVectorIterator(const GraphElementVectorIterator<!IsConst, IsReverse, T, Alloc> & other)
            : m_it(other.m_it)
            {
                static_assert(IsConst, "Can't assign non-const GraphElementVectorIterator from const GraphElementVectorIterator");
            }
            
            GraphElementVectorIterator &  operator=(const GraphElementVectorIterator & other)
            {
                if (this != &other) {
                    m_it = other.m_it;
                }
                return *this;
            }
            
            GraphElementVectorIterator & operator=(const GraphElementVectorIterator<!IsConst, IsReverse, T, Alloc> & other)
            {
                static_assert(IsConst, "Can't assign non-const GraphElementVectorIterator from const GraphElementVectorIterator");
                m_it = other.m_it;
                return *this;
            }
            
            reference operator*()
            {
                return *m_it;
            }
            
            pointer operator->()
            {
                return *m_it;
            }
            
            bool operator ==(const GraphElementVectorIterator & other) const
            {
                return m_it == other.m_it;
            }
            
            bool operator !=(const GraphElementVectorIterator & other) const
            {
                return m_it != other.m_it;
            }
            
            bool operator <(const GraphElementVectorIterator & other) const
            {
                return m_it < other.m_it;
            }
            
            bool operator >(const GraphElementVectorIterator & other) const
            {
                return m_it > other.m_it;
            }
            
            GraphElementVectorIterator & operator++()
            {
                ++m_it;
                return *this;
            }
            
            GraphElementVectorIterator & operator--()
            {
                --m_it;
                return *this;
            }
            
            GraphElementVectorIterator operator++(int i)
            {
                auto it = *this;
                ++m_it;
                return it;
            }
            
            GraphElementVectorIterator operator--(int i)
            {
                auto it = *this;
                --m_it;
                return it;
            }
            
            
            GraphElementVectorIterator & operator+=(difference_type i)
            {
                m_it += i;
                return *this;
            }
            
            GraphElementVectorIterator & operator-=(difference_type i)
            {
                m_it -= i;
                return *this;
            }
            
            GraphElementVectorIterator operator+(difference_type i)
            {
                auto it = *this;
                it += i;
                return it;
            }
            
            GraphElementVectorIterator operator-(difference_type i)
            {
                auto it = *this;
                it += i;
                return it;
            }
            
            difference_type operator-(GraphElementVectorIterator & other)
            {
                return m_it - other.m_it;
            }
            
            reference operator[](size_t i)
            {
                return m_it[i];
            }
            
        private:
            
            Iterator m_it;
        };

        struct GraphElementVectorAccessor;
    }

    template <typename T, typename Alloc>
    class GraphElementVector
    {
        typedef std::vector<T*, Alloc> Vector;
        typedef std::vector<T const*, Alloc> ConstVector;
        friend struct Detail::GraphElementVectorAccessor;

    public:
        
        typedef typename Vector::value_type                                                     value_type;     
        typedef typename Vector::allocator_type                                                 allocator_type;
        typedef typename Vector::const_reference                                                reference;     
        typedef typename std::remove_reference<typename ConstVector::const_reference>::type     const_reference;  
        typedef typename Vector::const_pointer                                                  pointer;
        typedef typename ConstVector::const_pointer                                             const_pointer;
        typedef Detail::GraphElementVectorIterator<false, false, T, Alloc>                      iterator;
        typedef Detail::GraphElementVectorIterator<true, false, T, Alloc>                       const_iterator;
        typedef Detail::GraphElementVectorIterator<false, true, T, Alloc>                       reverse_iterator;
        typedef Detail::GraphElementVectorIterator<true, true, T, Alloc>                        const_reverse_iterator;
        typedef typename Vector::difference_type                                                difference_type;
        typedef typename Vector::size_type                                                      size_type;     
               
        explicit GraphElementVector(const allocator_type & alloc = allocator_type())
        : vec(alloc)
        {
        }
        
        explicit GraphElementVector(size_type n)
        : vec(n, NULL)
        {
        }
        
        GraphElementVector(size_type n, const value_type& val, const allocator_type& alloc = allocator_type())
        : vec(n, val, alloc)
        {
        }
        
        template <class InputIterator>
        GraphElementVector (InputIterator first, InputIterator last, const allocator_type& alloc = allocator_type())
        : vec(first, last, alloc)
        {
        }
        
        GraphElementVector(const GraphElementVector& x)
        : vec(x.vec)
        {
        }
        
        GraphElementVector(const GraphElementVector& x, const allocator_type& alloc)
        : vec(x.vec, alloc)
        {
        }
        
        GraphElementVector(GraphElementVector&& x)
        : vec(std::move(x.vec))
        {
        }
        
        GraphElementVector(GraphElementVector&& x, const allocator_type& alloc)
        : vec(std::move(x.vec), alloc)
        {
        }
        
        GraphElementVector(std::initializer_list<value_type> il, const allocator_type& alloc = allocator_type())
        : vec(il, alloc)
        {
        }
        
        ~GraphElementVector()
        {
        }
        
        friend bool operator== (const Gulo::GraphElementVector<T,Alloc>& lhs, const Gulo::GraphElementVector<T,Alloc>& rhs)
        {
            return lhs.vec == rhs.vec;
        }
        
        friend bool operator!= (const Gulo::GraphElementVector<T,Alloc>& lhs, const Gulo::GraphElementVector<T,Alloc>& rhs)
        {
            return lhs.vec != rhs.vec;
        }
        
        friend bool operator<= (const Gulo::GraphElementVector<T,Alloc>& lhs, const Gulo::GraphElementVector<T,Alloc>& rhs)
        {
            return lhs.vec <= rhs.vec;
        }
        
        friend bool operator>= (const Gulo::GraphElementVector<T,Alloc>& lhs, const Gulo::GraphElementVector<T,Alloc>& rhs)
        {
            return lhs.vec >= rhs.vec;
        }
        
        friend bool operator< (const Gulo::GraphElementVector<T,Alloc>& lhs, const Gulo::GraphElementVector<T,Alloc>& rhs)
        {
            return lhs.vec < rhs.vec;
        }
        
        friend bool operator> (const Gulo::GraphElementVector<T,Alloc>& lhs, const Gulo::GraphElementVector<T,Alloc>& rhs)
        {
            return lhs.vec > rhs.vec;
        }
        
        iterator begin() noexcept
        {
            return iterator(vec.begin());
        }
        
        const_iterator begin() const noexcept
        {
            return cbegin();
        }
        
        iterator end() noexcept
        {
            return iterator(vec.end());
        }
        
        const_iterator end() const noexcept
        {
            return cend();
        }
        
        reverse_iterator rbegin() noexcept
        {
            return reverse_iterator(vec.rbegin());
        }
        
        const_reverse_iterator rbegin() const noexcept
        {
            return crbegin();
        }
        
        reverse_iterator rend() noexcept
        {
            return reverse_iterator(vec.rend());
        }
        
        const_reverse_iterator rend() const noexcept
        {
            return crend();
        }
        
        const_iterator cbegin() const noexcept
        {
            return const_iterator(vec.cbegin());
        }
        
        const_iterator cend() const noexcept
        {
            return const_iterator(vec.cend());
        }
        
        const_reverse_iterator crbegin() const noexcept
        {
            return const_reverse_iterator(vec.crbegin());
        }
        
        const_reverse_iterator crend() const noexcept
        {
            return const_reverse_iterator(vec.crend());
        }
        
        size_type size() const noexcept
        {
            return vec.size();
        }
        
        size_type max_size() const noexcept
        {
            return vec.max_size();
        }
        
        size_type capacity() const noexcept
        {
            return vec.capacity();
        }
        
        bool empty() const noexcept
        {
            return vec.empty();
        }
        
        void reserve (size_type n)
        {
            vec.reserve(n);
        }
        
        void shrink_to_fit()
        {
            vec.shrink_to_fit();
        }
        
        reference operator[] (size_type n)
        {
            return vec[n];
        }
        
        const_reference operator[] (size_type n) const
        {
            return vec[n];
        }
        
        reference at (size_type n)
        {
            return vec.at(n);
        }
        
        const_reference at (size_type n) const
        {
            return vec.at(n);
        }
        
        reference front()
        {
            return vec.front();
        }
        
        const_reference front() const
        {
            return vec.front();
        }
        
        reference back()
        {
            return vec.back();
        }
        
        const_reference back() const
        {
            return vec.back();
        }
        
        value_type* data() noexcept
        {
            return vec.data();
        }
        
        const T * const * data() const noexcept
        {
            return vec.data();
        }

    protected:
        
        std::vector<T*, Alloc> & vector()
        {
            return vec;
        }
        
        const std::vector<T*, Alloc> & vector() const
        {
            return vec;
        }
        
        std::vector<T*, Alloc> vec;
    };
}

#endif // GULO_GRAPH_ELEMENTVECTOR_HPP
