#ifndef GULO_SEQIO_GENBANK_H
#define GULO_SEQIO_GENBANK_H

#include <algorithm>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/erase.hpp>
#include <fstream>                                                                         
#include "tokenizedsequencesource.h"

namespace Gulo
{
    namespace Detail
    {
        class GenbankSequenceParser : public SequenceSourceIteratorData
        {
        public:

            GenbankSequenceParser() : SequenceSourceIteratorData(), data(0), len(0) {}

            static bool isBadSeqChar(char ch) { return isspace(ch) || isdigit(ch);}
        	
            static const char * token() {return "LOCUS";}
            
            void setSource(std::pair<const char*, int> source)
            {
                data = source.first;
                len = source.second;
            }

            std::string definition()
            {
                static const char * definition = "\nDEFINITION";
                static const size_t definition_size = strlen(definition);
                static const char * brk = "\n";
                static const size_t brk_size = strlen(brk);

                const char * p = std::search(data + 1, data + len, definition, definition + definition_size);
                const char * q = std::search(p + definition_size, data+len, brk, brk + brk_size);
                while (q != data+len && *(q+1) != ' ')
                    q = std::search(q+1, data+len, brk, brk + brk_size);
                std::string ba(p + definition_size, q-p);
                boost::trim(ba);
                boost::erase_all(ba, "\n");
                return ba;
            }

            std::string sequence()
            {
                static const char * origin = "\nORIGIN ";
                static const size_t origin_size = strlen(origin);
                static const char * seqend = "//";
                static const size_t seqend_size = strlen(seqend);

                const char * p = std::search(data, data + len, origin, origin + origin_size);
                const char * q = std::search(p + origin_size, data+len, seqend, seqend + seqend_size);

                std::string ba(p + origin_size, q - p - origin_size);
                std::string::iterator it = std::remove_if(ba.begin(), ba.end(), isBadSeqChar);
                return ba.substr(0, it - ba.begin());
            }

        protected:

            const char * endOfLabel;
            const char * data;
            int len;
        };
    }

    struct GenbankSequenceSource : public Detail::TokenizedSequenceSource<Detail::GenbankSequenceParser>
    {
        GenbankSequenceSource(const char * filePath) : Detail::TokenizedSequenceSource<Detail::GenbankSequenceParser>(filePath) 
        {
        }

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
                if (line.size() > 5 && line.substr(0,5) == "LOCUS") {
                    f.close();
                    return true;
                }
            }
            return false;
        }
    };
}

#endif // SEQIO_GENBANK_H