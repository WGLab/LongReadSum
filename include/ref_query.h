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
        uint32_t cpg_site_count = 0;

        // Map of reference position (0-indexed) to CpG site (true/false) on the
        // forward strand for each chromosome
        // std::map<std::string, std::map<int64_t, bool>> chr_pos_to_cpg;

        // Map of reference position (0-indexed) to CpG site (true/false) on the
        // reverse strand for each chromosome
        // std::map<std::string, std::map<int64_t, bool>> chr_pos_to_cpg_rev;

        // Map of reference position (0-indexed) to CpG site (true/false) for
        // all chromosomes
        std::map<std::string, std::map<int64_t, bool>> chr_pos_to_cpg;

        // Reverse strand CpG site map
        // std::map<std::string, std::map<int64_t, bool>> chr_pos_to_cpg_rev;

    public:
        int setFilepath(std::string fasta_filepath);
        std::string getFilepath();
        char getBase(std::string chr, int64_t pos);
        void generateCpGMap();
        void addCpGSiteModification(std::string chr, int64_t pos, int strand);
        std::pair<uint64_t, uint64_t> getCpGModificationCounts(int strand);
};

#endif // REF_QUERY_H
