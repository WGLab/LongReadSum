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

        // Map of reference position (0-indexed) to CpG site (true/false)
        std::map<std::string, std::map<int64_t, bool>> chr_pos_to_cpg;

    public:
        int setFilepath(std::string fasta_filepath);
        std::string getFilepath();
        char getBase(std::string chr, int64_t pos);
        void generateCpGMap();
        bool isCpG(std::string chr, int64_t pos);
};

#endif // REF_QUERY_H
