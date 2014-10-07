#ifndef GULO_SEQIO_TOKENIZEDSEQUENCESOURCE_H
#define GULO_SEQIO_TOKENIZEDSEQUENCESOURCE_H

// for MSVC, may need to patch boost: see https://svn.boost.org/trac/boost/ticket/9332

#include "sequencesource.h"
#include <boost/cstdint.hpp>
#include <boost/interprocess/file_mapping.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <cstring>

namespace Gulo
{
    namespace Detail
    {
        class TokenizedSequenceSourceView : public SequenceSourceView
        {
        public:

            enum file_pos
            {
                begin,
                end
            };

            TokenizedSequenceSourceView(boost::interprocess::file_mapping * fileMapping, const char * token, file_pos pos = begin, boost::intmax_t page_size = 64000000);
            TokenizedSequenceSourceView(const TokenizedSequenceSourceView & other);
            ~TokenizedSequenceSourceView();

            TokenizedSequenceSourceView & operator =(const TokenizedSequenceSourceView & other);

            bool operator ==(const SequenceSourceView& other) const;
            bool operator !=(const SequenceSourceView& other) const {return !(*this == other);}

             bool atEnd() const {return pos1 == file_size;}
             bool atRend() const {return pos1 == -1;}
             std::pair<const char*, int> next();
             std::pair<const char*, int> prev();

        private:

            boost::interprocess::file_mapping * file;
            boost::intmax_t file_size;
            boost::intmax_t page_size;
            boost::intmax_t cur_page_size;
            boost::intmax_t offset;
            boost::intmax_t pos1;
            boost::intmax_t pos2;
        	boost::interprocess::mapped_region * region;
            const char * buffer;
            const char * token;
			boost::intmax_t token_size;
        };

        template <typename Parser>
        class TokenizedSequenceSource : public SequenceSource
        {
        public:

        	 TokenizedSequenceSource(const char * filePath) : SequenceSource(filePath) {map = new boost::interprocess::file_mapping(filePath, boost::interprocess::read_only);}
        	 ~TokenizedSequenceSource() {delete map;}

             iterator begin() {iterator it(new TokenizedSequenceSourceView(map, Parser::token(), TokenizedSequenceSourceView::begin, 64000000), new Parser()); ++it; return it;}
             iterator end() {return iterator(new TokenizedSequenceSourceView(map, Parser::token(), TokenizedSequenceSourceView::end, 1), new Parser());}
             reverse_iterator rbegin() {reverse_iterator it(new TokenizedSequenceSourceView(map, Parser::token(), TokenizedSequenceSourceView::end, 64000000), new Parser()); ++it; return it;}
             reverse_iterator rend() {return reverse_iterator(new TokenizedSequenceSourceView(map, Parser::token(), TokenizedSequenceSourceView::begin, 1), new Parser());}

        private:

        	boost::interprocess::file_mapping * map;
        };

        TokenizedSequenceSourceView::TokenizedSequenceSourceView(boost::interprocess::file_mapping * _file, const char * _token, file_pos pos, boost::intmax_t _page_size) :
            file(_file),
            file_size(boost::filesystem::file_size(file->get_name())),
            page_size(_page_size),
            cur_page_size(std::min(page_size, file_size)),
            offset((pos == begin) ? 0 : file_size - cur_page_size),
            pos1((pos == begin) ? -1 : file_size),
            pos2((pos == begin) ? -1 : file_size),
            region(new boost::interprocess::mapped_region(*file, file->get_mode(), offset, cur_page_size)),
            buffer((const char*)(region->get_address())),
            token(_token),
            token_size(strlen(token))
        {
        }

        TokenizedSequenceSourceView::TokenizedSequenceSourceView(const TokenizedSequenceSourceView & other) :
            file(other.file),
            file_size(other.file_size),
            page_size(other.page_size),
            cur_page_size(other.cur_page_size),
            offset(other.offset),
            pos1(other.pos1),
            pos2(other.pos2),
            region(new boost::interprocess::mapped_region(*file, file->get_mode(), offset, cur_page_size)),
            buffer((const char*)(region->get_address())),
            token(other.token),
            token_size(other.token_size)
        {
        }

        TokenizedSequenceSourceView::~TokenizedSequenceSourceView()
        {
            delete region;
        }

        TokenizedSequenceSourceView & TokenizedSequenceSourceView::operator =(const TokenizedSequenceSourceView & other)
        {
            if (this != &other) {
                file = other.file;
                page_size = other.page_size;
                cur_page_size = other.cur_page_size;
                file_size = other.file_size;
                pos1 = other.pos1;
                pos2 = other.pos2;
                token = other.token;
                token_size = other.token_size;
                offset = other.offset;
                delete region;
                region = new boost::interprocess::mapped_region(*file, file->get_mode(), offset, cur_page_size);
                buffer = (const char*)(region->get_address());
            }
            return *this;
        }

        bool TokenizedSequenceSourceView::operator ==(const SequenceSourceView & o) const
        {
            TokenizedSequenceSourceView & other = (TokenizedSequenceSourceView&)o;
            return ((file == other.file) && (pos1 == other.pos1) && (pos2 == other.pos2)/* && strcmp(token, other.token) == 0*/);
        }

        std::pair<const char*, int> TokenizedSequenceSourceView::next()
        {
            if (pos1 < file_size || pos1 == -1)
                pos1 = pos2 + 1;
            else {
                return std::pair<const char*, int>((const char*) 0, 0);
            }

            do {
                pos1 = std::search(buffer + pos1 - offset, buffer + cur_page_size, token, token + token_size) - buffer + offset;
                while (pos1 < offset + cur_page_size && !((pos1 == 0 || *(buffer + pos1 - 1 - offset) == '\n'))) {
                    ++pos1;
                    pos1 = std::search(buffer + pos1 - offset, buffer + cur_page_size, token, token + token_size) - buffer + offset;
                }
                if (pos1 >= cur_page_size + offset && pos1 < file_size) {
                    cur_page_size = std::min(file_size - (pos1 - token_size + 1), page_size);
                    offset = pos1 = pos1 - token_size + 1;
                    delete region;
                    region = new boost::interprocess::mapped_region(*file, file->get_mode(), offset, cur_page_size);
                    buffer = (const char*)(region->get_address());
                } else {
                    break;
                }
            } while (true);

            if (pos1 == file_size) {
                pos2 = pos1;
                return std::pair<const char*, int>((const char*)0, 0);
            }

            pos2 = pos1 + 1;

            do {
                pos2 = std::search(buffer + pos2 - offset, buffer + cur_page_size, token, token + token_size) - buffer + offset;
                while (pos2 < offset + cur_page_size && !((pos2 == 0 || *(buffer + pos2 - 1 - offset) == '\n'))) {
                    ++pos2;
                    pos2 = std::search(buffer + pos2 - offset, buffer + cur_page_size, token, token + token_size) - buffer + offset;
                }
                if (pos2 == cur_page_size + offset && pos2 < file_size) {
                    cur_page_size = std::min(file_size - (pos2 - token_size + 1), page_size);
                    offset = pos2 = pos2 - token_size + 1;
                    delete region;
                    region = new boost::interprocess::mapped_region(*file, file->get_mode(), offset, cur_page_size);
                    buffer = (const char*)(region->get_address());
                } else {
                    --pos2;
                    break;
                }
            } while (true);

            if (pos1 < offset) {
                cur_page_size = std::max(file_size - pos1, std::max(page_size, pos2 - pos1));
                offset = pos1;
                delete region;
                region = new boost::interprocess::mapped_region(*file, file->get_mode(), offset, cur_page_size);
                buffer = (const char*)(region->get_address());
            }

            return std::pair<const char*, int>(buffer - offset + pos1, (int)(pos2 -  pos1));
        }

        std::pair<const char*, int> TokenizedSequenceSourceView::prev()
        {
            if (pos2 != -1 && pos1 != -1 && pos1 != 0)
                pos2 = pos1 - 1;
            else {
                pos2 = pos1 = -1;
                return std::pair<const char*, int>((const char*) 0, 0);
            }

            if (pos2 < offset) {
                cur_page_size = std::min(pos2, page_size);
                offset = pos2 - cur_page_size;
                delete region;
                region = new boost::interprocess::mapped_region(*file, file->get_mode(), offset, cur_page_size);
                buffer = (const char*)(region->get_address());
            }

            pos1 = pos2 - 1;

            do {
                const char * max = buffer + pos1 - offset;
                const char * ch = std::find_end(buffer, max, token, token + token_size);
                pos1 = ch - buffer + offset;
                while (ch < max && !((pos1 == 0 || *(buffer + pos1 - 1 - offset) == '\n'))) {
                    --pos1;
                    max = buffer + pos1 - offset;
                    ch = std::find_end(buffer, max, token, token + token_size);
                    pos1 = ch - buffer + offset;
                }
                if ((ch == max) && pos1 > 0) {
                    cur_page_size = std::min(file_size, std::min(offset + token_size + 1, page_size + token_size + 1));
                    offset = (cur_page_size > offset) ? 0 : offset - cur_page_size + 1;
                    pos1 = offset + cur_page_size;
                    delete region;
                    region = new boost::interprocess::mapped_region(*file, file->get_mode(), offset, cur_page_size);
                    buffer = (const char*)(region->get_address());
                } else {
                    break;
                }
            } while (true);

            if (pos1 == -1) {
                pos2 = pos1 = -1;
                return std::pair<const char*, int>((const char*) 0, 0);
            }

            if (pos2 >= offset + cur_page_size) {
                cur_page_size = std::min(std::max(page_size, (pos2 - offset) + 1), file_size - pos2 + 1);
                offset = (cur_page_size > pos2) ? 0 : pos2 - cur_page_size;
                delete region;
                region = new boost::interprocess::mapped_region(*file, file->get_mode(), offset, cur_page_size);
                buffer = (const char*)(region->get_address());
            }

            return std::pair<const char*, int>(buffer + pos1 - offset, (int)(pos2 -  pos1));
        }
    }
}

#endif // GULO_SEQIO_TOKENIZEDSEQUENCESOURCE_H
