#ifndef GULO_SEQIO_PHYLIP_H
#define GULO_SEQIO_PHYLIP_H

#include "sequencesource.h"
#include <vector>
#include <fstream>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/erase.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>

namespace Gulo
{
    namespace Detail
    {
        class PhylipSequenceParser : public SequenceSourceIteratorData
        {
        public:

            PhylipSequenceParser(const std::vector<std::pair<std::string, std::string> > & vector) : SequenceSourceIteratorData(), pos(0), vec(vector) {}

        	void setSource(std::pair<const char*, int> data);
            std::string definition();
            std::string sequence();

        protected:

        	int pos;
        	const std::vector<std::pair<std::string, std::string> > & vec;
        };

        struct PhylipSequenceSourceView : public SequenceSourceView
        {
            PhylipSequenceSourceView(size_t iter, const std::vector<std::pair<std::string, std::string> > & vector) : it(iter), vec(vector) {}

            bool atEnd() const;
            bool atRend() const;
            std::pair<const char*, int> next();
            std::pair<const char*, int> prev();
            bool operator ==(const SequenceSourceView&) const;
            bool operator !=(const SequenceSourceView&) const;

        protected:

        	size_t it;
        	const std::vector<std::pair<std::string, std::string> > & vec;
        };
    }

    struct PhylipSequenceSource : public SequenceSource
    {
        PhylipSequenceSource(const char * filePath);

    	static bool isValidFile(const char * filePath);

    	iterator begin();
        iterator end();
        reverse_iterator rbegin();
        reverse_iterator rend();

    protected:

        int meta_seqLen;
        int meta_seqCount;
    	std::vector<std::pair<std::string,std::string> > data;
    };

    bool PhylipSequenceSource::isValidFile(const char * filePath)
    {
        //bool gotLabels = false;
        std::string line;
        std::ifstream myfile (filePath);
        if (myfile.is_open())
        {
            while (myfile.good() ) {
                getline (myfile,line);
                boost::trim(line);
                if (line.size() > 0) {
                    std::vector<std::string> splitVec;
                    boost::split(splitVec, line, boost::is_space(), boost::token_compress_on );
                    if (splitVec.size() == 2) {
                        try
                        {
                            boost::lexical_cast<int>(splitVec[0]);
                        }
                        catch(boost::bad_lexical_cast &)
                        {
                            myfile.close();
                            return false;
                        }
                        try
                        {
                            boost::lexical_cast<int>(splitVec[1]);
                        }
                        catch(boost::bad_lexical_cast &)
                        {
                            myfile.close();
                            return false;
                        }
                        myfile.close();
                        return true;
                    } else {
                        myfile.close();
                        return false;
                    }
                } 
            }
            myfile.close();
            return false;
        }
        return false;
    }

    PhylipSequenceSource::PhylipSequenceSource(const char * filePath) : 
        SequenceSource(filePath) ,
        meta_seqLen(-1),
        meta_seqCount(-1)
    {
        bool gotLabels = false;
        std::string line;
        std::ifstream myfile (filePath);
        meta_seqCount = -1;
        meta_seqLen = -1;
        unsigned int pos = 0;
        if (myfile.is_open())
        {
            while ( myfile.good() ) {
                getline (myfile,line);
                boost::trim(line);
                if (meta_seqCount == -1 && line.size() > 0) {
                    std::vector<std::string> splitVec;
                    boost::split(splitVec, line, boost::is_space(), boost::token_compress_on );
                    if (splitVec.size() == 2) {
                        try
                        {
                            meta_seqCount = boost::lexical_cast<int>(splitVec[0]);
                        }
                        catch(boost::bad_lexical_cast &)
                        {
                            meta_seqCount = -1;
                            meta_seqLen = -1;
                        }
                        try
                        {
                            meta_seqLen = boost::lexical_cast<int>(splitVec[1]);
                        }
                        catch(boost::bad_lexical_cast &)
                        {
                            meta_seqCount = -1;
                            meta_seqLen = -1;
                        }
                    }
                } else if (meta_seqCount != -1) {
                    if (!gotLabels) {
                        data.push_back(std::pair<std::string,std::string>(line.substr(0, 10), line.substr(10, line.size()-10)));
                        ++pos;
                        if (pos >= (unsigned int) meta_seqCount) {
                            pos = 0;
                            gotLabels = true;
                        }
                    } else  {
                        data[pos].second += line;
                        ++pos;
                        if (pos >= data.size())
                            pos = 0;
                    }
                }
            }
            myfile.close();
        }

        for (unsigned int i = 0; i< data.size(); ++i) {
            boost::trim(data[i].first);
            boost::trim(data[i].second);
            boost::erase_all(data[i].second, " ");
        }
    }

    SequenceSource::iterator PhylipSequenceSource::begin()
    {
        SequenceSource::iterator it(new Detail::PhylipSequenceSourceView(-1, data), new Detail::PhylipSequenceParser(data));
        ++it;
        return it;
    }

    SequenceSource::iterator PhylipSequenceSource::end()
    {
        return SequenceSource::iterator(new Detail::PhylipSequenceSourceView(data.size(), data), new Detail::PhylipSequenceParser(data));
    }

    SequenceSource::reverse_iterator PhylipSequenceSource::rbegin()
    {
        SequenceSource::reverse_iterator it(new Detail::PhylipSequenceSourceView(data.size() - 1, data), new Detail::PhylipSequenceParser(data));
        ++it;
        return it;
    }

    SequenceSource::reverse_iterator PhylipSequenceSource::rend()
    {
        return SequenceSource::reverse_iterator(new Detail::PhylipSequenceSourceView(-1, data), new Detail::PhylipSequenceParser(data));
    }

    namespace Detail
    {

        bool PhylipSequenceSourceView::atEnd() const
        {
            return (it == vec.size());
        }

        bool PhylipSequenceSourceView::atRend() const
        {
            return (it == static_cast<size_t>(-1));
        }

        std::pair<const char*, int> PhylipSequenceSourceView::next()
        {
            ++it;
            return std::pair<const char*, int>((const char*)0,it);
        }

        std::pair<const char*, int> PhylipSequenceSourceView::prev()
        {
            --it;
            return std::pair<const char*, int>((const char*)0,it);
        }

        bool PhylipSequenceSourceView::operator ==(const SequenceSourceView& other) const
        {
            return it == ((const PhylipSequenceSourceView&)(other)).it;
        }

        bool PhylipSequenceSourceView::operator !=(const SequenceSourceView& other) const
        {
            return it != ((const PhylipSequenceSourceView&)(other)).it;
        }

        void PhylipSequenceParser::setSource(std::pair<const char*, int> data)
        {
            pos = data.second;
        }

        std::string PhylipSequenceParser::definition()
        {
            return vec[pos].first;
        }

        std::string PhylipSequenceParser::sequence()
        {
            return vec[pos].second;
        }
    }

}

#endif // GULO_SEQIO_PHYLIP_H