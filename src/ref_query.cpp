#include "ref_query.h"

#include <string.h>
#include <iostream>
#include <unordered_map>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>


int RefQuery::setFilepath(std::string fasta_filepath)
{
    if (fasta_filepath == "")
    {
        std::cout << "No FASTA filepath provided" << std::endl;
        return 1;
    }

    this->fasta_filepath = fasta_filepath;

    // Parse the FASTA file
    std::ifstream fasta_file(fasta_filepath);
    if (!fasta_file.is_open())
    {
        std::cout << "Could not open FASTA file " << fasta_filepath << std::endl;
        exit(1);
    }

    // Get the chromosomes and sequences
    std::vector<std::string> chromosomes;
    std::unordered_map<std::string, std::string> chr_to_seq;
    std::string current_chr = "";
    std::string sequence = "";
    std::string line_str = "";
    while (std::getline(fasta_file, line_str))
    {
        // Check if the line is a header
        if (line_str[0] == '>')
        {
            // Header line, indicating a new chromosome
            // Store the previous chromosome and sequence
            if (current_chr != "")
            {
                chromosomes.push_back(current_chr);  // Add the chromosome to the list
                chr_to_seq[current_chr] = sequence;  // Add the sequence to the map
                sequence = "";  // Reset the sequence
            }

            // Get the new chromosome
            current_chr = line_str.substr(1);

            // Remove the description
            size_t space_pos = current_chr.find(" ");
            if (space_pos != std::string::npos)
            {
                current_chr.erase(space_pos);
            }

            // Check if the chromosome is already in the map
            if (chr_to_seq.find(current_chr) != chr_to_seq.end())
            {
                std::cerr << "Duplicate chromosome " << current_chr << std::endl;
                exit(1);
            }
        } else {
            // Sequence line
            sequence += line_str;
        }
    }

    // Add the last chromosome at the end of the file
    if (current_chr != "")
    {
        chromosomes.push_back(current_chr);  // Add the chromosome to the list
        chr_to_seq[current_chr] = sequence;  // Add the sequence to the map
    }

    // Close the file
    fasta_file.close();

    // Sort the chromosomes
    std::sort(chromosomes.begin(), chromosomes.end());

    // Set the chromosomes and sequences
    this->chromosomes = chromosomes;
    this->chr_to_seq = chr_to_seq;

    // Find CpG sites
    this->generateCpGMap();

    return 0;
}

void RefQuery::generateCpGMap()
{    
    // Iterate over each chromosome
    std::cout << "Locating CpG sites..." << std::endl;
    for (const std::string& chr : this->chromosomes)
    {
        this->chr_to_cpg[chr] = std::unordered_set<int64_t>();
        this->chr_to_cpg_mod[chr] = std::unordered_set<int64_t>();

        // Iterate over each position in the sequence
        const std::string& sequence = this->chr_to_seq[chr];
        for (int32_t pos = 0; pos < (int32_t)sequence.size(); pos++)
        {
            // Check if the base is a C
            if (sequence[pos] == 'C')
            {
                // Check if the next base is a G
                if (pos + 1 < (int32_t)sequence.size() && sequence[pos + 1] == 'G')
                {
                    // Add the CpG site to the map (1-based index)
                    int32_t pos1 = pos + 1;
//                    this->chr_pos_to_cpg[chr][pos1] = false;  // Initialize as false since no modifications have been found
                    this->cpg_total_count++;

                    this->chr_to_cpg[chr].insert(pos1);

                    // Skip the next base
                    // pos++;
                }
            }
        }
    }
}

void RefQuery::addCpGSiteModification(std::string chr, int64_t pos, int strand)
{
    // Update the CpG site if it exists
    // Reverse strand (position is the G in the CpG site, so move back one
    // position to get the C position stored in the map)
    if (strand == 1) {
        pos--;
    }

    // Find the CpG site in the map
//    if (this->chr_pos_to_cpg[chr].find(pos) != this->chr_pos_to_cpg[chr].end())
    if (this->chr_to_cpg[chr].find(pos) != this->chr_to_cpg[chr].end())
    {
        this->test_count++;

        // Update the CpG site if not already modified
        if (this->chr_to_cpg_mod[chr].find(pos) == this->chr_to_cpg_mod[chr].end())
        {
            this->chr_to_cpg_mod[chr].insert(pos);
            this->cpg_modified_count++;
        }

        // Update the modified count
//        this->cpg_modified_count++;

        // Update the CpG site if not already modified
//        if (!this->chr_pos_to_cpg[chr][pos])
//        {
//            this->chr_pos_to_cpg[chr][pos] = true;
//            this->cpg_modified_count++;
//        }
//        this->chr_pos_to_cpg[chr][pos] = true;
//        this->cpg_modified_count++;
    }
}

std::pair<uint64_t, uint64_t> RefQuery::getCpGModificationCounts(int strand)
{
    uint64_t modified_count = 0;
    uint64_t unmodified_count = 0;

    std::cout << "Calculating CpG modification counts..." << std::endl;
    modified_count = this->cpg_modified_count;
    unmodified_count = this->cpg_total_count - this->cpg_modified_count;
    std::cout << "Modified CpG sites: " << modified_count << std::endl;
    std::cout << "Unmodified CpG sites: " << unmodified_count << std::endl;
    std::cout << "Total CpG sites: " << modified_count + unmodified_count << std::endl;
    std::cout << "Percentage of CpG sites modified: " << (double)modified_count / (modified_count + unmodified_count) * 100 << "%" << std::endl;
    return std::make_pair(modified_count, unmodified_count);
}

std::string RefQuery::getFilepath()
{
    return this->fasta_filepath;
}

// Function to get the reference sequence at a given position range
std::string RefQuery::getBase(std::string chr, int64_t pos)
{
    // Convert positions from 1-indexed (reference) to 0-indexed (string indexing)
    pos--;

    // Ensure that the position is not negative
    if (pos < 0)
    {
        return "N";
    }

    // Get the sequence
    const std::string& sequence = this->chr_to_seq[chr];

    // Check if the position is out of range
    if (static_cast<std::size_t>(pos) >= sequence.size()) {
        return "N";
    }

    // Get the base
    std::string base = sequence.substr(static_cast<std::size_t>(pos), 1);
//    return str[static_cast<std::size_t>(index)];

    // Get the base
//    char base = sequence[pos];

    // If the base is empty, return empty string
    if (base == "")
    {
        return "N";
    }
//    if (base == '\0')
//    {
//        return 'N';
//    }

    return base;
}
