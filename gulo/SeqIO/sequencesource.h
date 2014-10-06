#ifndef GULO_SEQIO_SEQUENCESOURCE_H
#define GULO_SEQIO_SEQUENCESOURCE_H

#include <string>

namespace Gulo
{
    namespace Detail
    {
        struct SequenceSourceIteratorData
        {
            virtual ~SequenceSourceIteratorData() {}
            virtual std::string definition() {return "";}
            virtual std::string sequence() {return "";}
            virtual void setSource(std::pair<const char*, int> data) = 0;
        };

        struct SequenceSourceView
        {
            virtual ~SequenceSourceView() {}
            virtual bool atEnd() const = 0;
            virtual bool atRend() const = 0;
            virtual std::pair<const char*, int> next() = 0;
            virtual std::pair<const char*, int> prev() = 0;
            virtual bool operator ==(const SequenceSourceView&) const = 0;
            virtual bool operator !=(const SequenceSourceView&) const = 0;
        };

        template <bool IsReverse>
        struct SequenceSourceIterator
        {
            SequenceSourceIterator(SequenceSourceView * _view, SequenceSourceIteratorData * _data) : view(_view), data(_data) {}
            virtual ~SequenceSourceIterator() {delete data; delete view;}

            SequenceSourceIteratorData & operator *() {return *data;}
            SequenceSourceIteratorData * operator ->() {return data;}
            SequenceSourceIterator & operator++() {if (!IsReverse) data->setSource(view->next()); else data->setSource(view->prev()); return *this;}
            SequenceSourceIterator & operator--() {if (!IsReverse) data->setSource(view->prev()); else data->setSource(view->next()); return *this;}
            bool operator != (const SequenceSourceIterator & other) const {return (*view != *(other.view));}
            bool operator == (const SequenceSourceIterator & other) const {return (*view == *(other.view));}
        	bool atEnd() const {return IsReverse ? view->atRend() : view->atEnd();}

        protected:

            SequenceSourceView * view;
            SequenceSourceIteratorData * data;
        };
    } 

    struct SequenceSource
    {
        typedef Detail::SequenceSourceIterator<false> iterator;
        typedef Detail::SequenceSourceIterator<true> reverse_iterator;

        SequenceSource(const char * filePath) : path(filePath) {}
        virtual ~SequenceSource() {}

        virtual iterator begin() = 0;
        virtual iterator end() = 0;
        virtual reverse_iterator rbegin() = 0;
        virtual reverse_iterator rend() = 0;

    protected:

        const char * path;
    };
}

#endif // GULO_SEQIO_SEQUENCESOURCE_H
