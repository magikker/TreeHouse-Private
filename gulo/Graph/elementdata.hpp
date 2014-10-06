#ifndef GULO_GRAPH_ELEMENTDATA_HPP
#define GULO_GRAPH_ELEMENTDATA_HPP

namespace Gulo
{
    namespace Detail
    {
        template <typename DataType, typename CustomDataType>
        class GraphElementData
        {
        public:

            GraphElementData(const DataType & data = DataType(), const CustomDataType & customData = CustomDataType())
            : m_data(data),
              m_customData(customData)
            {
            }

            GraphElementData(DataType && data, CustomDataType && customData)
            : m_data(std::move(data)),
              m_customData(std::move(customData))
            {
            }

            GraphElementData(const GraphElementData & other)
            : m_data(other.m_data),
              m_customData(other.m_customData)
            {
            }

            GraphElementData(GraphElementData && other)
            : m_data(std::move(other.m_data)),
              m_customData(std::move(other.m_customData))
            {
            }

            virtual ~GraphElementData()
            {
            }

            GraphElementData & operator=(const GraphElementData & other)
            {
                if (this != &other) {
                    m_data = other.m_data;
                    m_customData = other.m_customData;
                }
                return *this;
            }

            GraphElementData & operator=(GraphElementData && other)
            {
                m_data = std::move(other.m_data);
                m_customData = std::move(other.m_customData);
                return *this;
            }

            DataType & value()
            {
                return m_data;
            }
            
            const DataType & value() const
            {
                return m_data;
            }

            CustomDataType & customValue()
            {
                return m_customData;
            }

            const CustomDataType & customValue() const
            {
                return m_customData;
            }

            constexpr static bool storesData()
            {
                return true;
            }

            constexpr static bool storesCustomData()
            {
                return true;
            }

        protected:

            DataType m_data;
            CustomDataType m_customData;
        };

        template <typename DataType>
        class GraphElementData<DataType, void>
        {
        public:

            GraphElementData(const DataType & data = DataType())
            : m_data(data)
            {
            }

            GraphElementData(DataType && data)
            : m_data(std::move(data))
            {
            }

            GraphElementData(const GraphElementData & other)
            : m_data(other.m_data)
            {
            }

            GraphElementData(GraphElementData && other)
            : m_data(std::move(other.m_data))
            {
            }

            virtual ~GraphElementData()
            {
            }

            GraphElementData & operator=(const GraphElementData & other)
            {
                if (this != &other) {
                    m_data = other.m_data;
                }
                return *this;
            }

            GraphElementData & operator=(GraphElementData && other)
            {
                m_data = std::move(other.m_data);
                return *this;
            }

            DataType & value()
            {
                return m_data;
            }
            
            const DataType & value() const
            {
                return m_data;
            }

            constexpr static bool storesData()
            {
                return true;
            }

            constexpr static bool storesCustomData()
            {
                return false;
            }

        protected:

            DataType m_data;
        };

        template <typename CustomDataType>
        class GraphElementData<void, CustomDataType>
        {
        public:

            GraphElementData(const CustomDataType & customData = CustomDataType())
            : m_customData(customData)
            {
            }

            GraphElementData(CustomDataType && customData)
            : m_customData(std::move(customData))
            {
            }

            GraphElementData(const GraphElementData & other)
            : m_customData(other.m_customData)
            {
            }

            GraphElementData(GraphElementData && other)
            : m_customData(std::move(other.m_customData))
            {
            }

            virtual ~GraphElementData()
            {
            }

            GraphElementData & operator=(const GraphElementData & other)
            {
                if (this != &other) {
                    m_customData = other.m_customData;
                }
                return *this;
            }

            GraphElementData & operator=(GraphElementData && other)
            {
                m_customData = std::move(other.m_customData);
                return *this;
            }

            CustomDataType & customValue()
            {
                return m_customData;
            }

            const CustomDataType & customValue() const
            {
                return m_customData;
            }

            constexpr static bool storesData()
            {
                return false;
            }

            constexpr static bool storesCustomData()
            {
                return true;
            }

        protected:

            CustomDataType m_customData;
        };

        template <>
        class GraphElementData<void, void>
        {
        public:

            constexpr static bool storesData()
            {
                return false;
            }

            constexpr static bool storesCustomData()
            {
                return false;
            }
        };
    }
}

#endif // GULO_GRAPH_ELEMENTDATA_HPP
