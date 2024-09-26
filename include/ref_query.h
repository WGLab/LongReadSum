// RefQuery: A class for querying a reference genome in FASTA format

#ifndef REF_QUERY_H
#define REF_QUERY_H

#include <stdint.h>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>

class RefQuery {
    private:
        std::string fasta_filepath;
        std::vector<std::string> chromosomes;
        std::unordered_map<std::string, std::string> chr_to_seq;
        uint64_t cpg_modified_count = 0;
        uint64_t cpg_total_count = 0;

        // Map of chromosome to CpG site positions
        std::unordered_map<std::string, std::unordered_set<int64_t>> chr_to_cpg;

        // Map of chromosome to CpG site positions with modifications
        std::unordered_map<std::string, std::unordered_set<int64_t>> chr_to_cpg_mod;

    public:
        int setFilepath(std::string fasta_filepath);
        std::string getFilepath();
        std::string getBase(std::string chr, int64_t pos);
        void generateCpGMap();
        void addCpGSiteModification(std::string chr, int64_t pos, int strand);
        std::pair<uint64_t, uint64_t> getCpGModificationCounts(int strand);
};

#endif // REF_QUERY_H
