#ifndef GULO_SEQIO_FASTA_H
#define GULO_SEQIO_FASTA_H
                                                                               
#include "tokenizedsequencesource.h"
#include <algorithm>
#include <boost/algorithm/string/trim.hpp>
#include <fstream>
#include <iostream>

namespace Gulo
{
    namespace Detail
    {
        class FastaSequenceParser : public SequenceSourceIteratorData
        {
        public:

            FastaSequenceParser() : SequenceSourceIteratorData(), endOfLabel(0), data(0), len(0) {}

        	static const char * token() 
            {
                return ">";
            }
            
            void setSource(std::pair<const char*, int> source)
            {
                endOfLabel = 0;
                data = source.first;
                len = source.second;
            }

            std::string definition()
            {
                if (endOfLabel == 0)
                    endOfLabel = std::find(data + 1, data + len, '\n');
                std::string result(data + 1, endOfLabel - data);
                boost::trim(result);
                return result;
            }

            std::string sequence()
            {
                if (endOfLabel == 0)
                    endOfLabel = std::find(data + 1, data + len, '\n');
                std::string result(endOfLabel + 1, len - (endOfLabel - data));
                result.erase( std::remove_if(result.begin(), result.end(), isspace), result.end() );
                boost::trim(result);
                return result;
            }

        protected:

            const char * endOfLabel;
            const char * data;
            int len;
        };
    }

    struct FastaSequenceSource : public Detail::TokenizedSequenceSource<Detail::FastaSequenceParser>
    {
        FastaSequenceSource(const char * filePath) : Detail::TokenizedSequenceSource<Detail::FastaSequenceParser>(filePath) {}
    	
        static bool isValidFile(const char * filePath)
        {
            std::string line;
            std::ifstream f(filePath);
            if (f.is_open()) {
                std::getline(f, line);
                boost::trim(line);
                while (line.length() == 0 && f.good()) {
                    std::getline(f, line);
                    boost::trim(line);
                }
                if (line.size() > 0 && line[0] == '>') {
                    f.close();
                    return true;
                }
            }
            return false;
        }

    };
}

#endif // GULO_SEQIO_FASTA_H