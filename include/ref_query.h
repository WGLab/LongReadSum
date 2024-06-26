// RefQuery: A class for querying a reference genome in FASTA format

#ifndef REF_QUERY_H
#define REF_QUERY_H

#include <string>
#include <map>
#include <unordered_map>
#include <vector>

class RefQuery {
    private:
        std::string fasta_filepath;
        std::vector<std::string> chromosomes;
        std::unordered_map<std::string, std::string> chr_to_seq;

    public:
        int setFilepath(std::string fasta_filepath);
        std::string getFilepath();
        char getBase(std::string chr, int64_t pos);
};

#endif // REF_QUERY_H
