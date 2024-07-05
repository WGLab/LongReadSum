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
        // Get the sequence
        const std::string& sequence = this->chr_to_seq[chr];

        // Iterate over each position in the sequence
        for (int64_t pos = 0; pos < (int64_t)sequence.size(); pos++)
        {
            // Check if the base is a C
            if (sequence[pos] == 'C')
            {
                // Check if the next base is a G
                if (pos + 1 < (int64_t)sequence.size() && sequence[pos + 1] == 'G')
                {
                    // Add the CpG site to the map
                    this->chr_pos_to_cpg[chr][pos] = false;  // Initialize as false since no modifications have been found

                    // Skip the next base
                    // pos++;

                    // Update the CpG site count
                    // this->cpg_site_count++;
                }
            }

            // // Check if the base is a G
            // else if (sequence[pos] == 'G')
            // {
            //     // Check if the previous base is a C
            //     if (pos - 1 >= 0 && sequence[pos - 1] == 'C')
            //     {
            //         // Add the CpG site to the map
            //         this->chr_pos_to_cpg_rev[chr][pos] = false;  // Initialize as false since no modifications have been found

            //         // // Update the CpG site count
            //         // this->cpg_site_count++;
            //     }
            // }
        }
    }
    std::cout << "CpG sites located." << std::endl;
}

void RefQuery::addCpGSiteModification(std::string chr, int64_t pos, int strand)
{
    // Update the CpG site
    this->chr_pos_to_cpg[chr][pos] = true;
    // Determine the map to use based on the strand
    // std::map<std::string, std::map<int64_t, bool>>* map;

    // // If on the forward strand, check if the position is a CpG site
    // if (strand == 0)
    // {
    //     map = &this->chr_pos_to_cpg;
    // } else {
    //     map = &this->chr_pos_to_cpg_rev;
    // }

    // // Update the CpG site
    // (*map)[chr][pos] = true;
}

std::pair<uint32_t, uint32_t> RefQuery::getCpGModificationCounts(int strand)
{
    uint32_t modified_count = 0;
    uint32_t unmodified_count = 0;

    std::cout << "Calculating CpG modification counts for strand " << strand << "..." << std::endl;

    // // Determine the map to use based on the strand
    // std::map<std::string, std::map<int64_t, bool>>* map;
    // if (strand == 0)
    // {
    //     map = &this->chr_pos_to_cpg;
    // } else {
    //     map = &this->chr_pos_to_cpg_rev;
    // }
    // Get the forward and reverse strand maps
    // std::map<std::string, std::map<int64_t, bool>>* forward_map = &this->chr_pos_to_cpg;
    // std::map<std::string, std::map<int64_t, bool>>* reverse_map =
    // &this->chr_pos_to_cpg_rev;
    
    // Iterate over each chromosome in the CpG site map
    for (const auto& chr_pos_map : this->chr_pos_to_cpg)
    {
        // Get the chromosome
        const std::string& chr = chr_pos_map.first;

        // Iterate over each CpG site in the chromosome
        for (const auto& pos_to_cpg : chr_pos_map.second)
        {
            // Get the position and CpG site
            int64_t pos = pos_to_cpg.first;
            bool is_cpg = pos_to_cpg.second;

            // Check if the CpG site is modified
            if (is_cpg)
            {
                // Increment the modified count
                modified_count++;
            } else {
                // Increment the unmodified count
                unmodified_count++;
            }
        }
    }


    // // Iterate over each chromosome in the forward map
    // for (const auto& chr_pos_map : *forward_map)
    // {
    //     // Get the chromosome
    //     const std::string& chr = chr_pos_map.first;

    //     // Iterate over each CpG site in the chromosome
    //     for (const auto& pos_to_cpg : chr_pos_map.second)
    //     {
    //         // Get the position and CpG site for the forward strand
    //         int64_t pos = pos_to_cpg.first;
    //         bool is_cpg = pos_to_cpg.second;

    //         // Get the CpG site for the reverse strand
    //         bool is_cpg_rev = false;
    //         if (reverse_map->find(chr) != reverse_map->end() && (*reverse_map)[chr].find(pos) != (*reverse_map)[chr].end())
    //         {
    //             is_cpg_rev = (*reverse_map)[chr][pos];
    //         }

    //         // Check if the CpG site is modified in either the forward or
    //         // reverse strand
    //         if (is_cpg || is_cpg_rev)
    //         {
    //             // Increment the modified count
    //             modified_count++;
    //         } else {
    //             // Increment the unmodified count
    //             unmodified_count++;
    //         }
    //     }
    // }

    // // Iterate over each chromosome in the map
    // for (const auto& chr_pos_map : *map)
    // {
    //     // Get the chromosome
    //     const std::string& chr = chr_pos_map.first;

    //     // Iterate over each CpG site in the chromosome
    //     for (const auto& pos_to_cpg : chr_pos_map.second)
    //     {
    //         // Get the position and CpG site
    //         int64_t pos = pos_to_cpg.first;
    //         bool is_cpg = pos_to_cpg.second;

    //         // Check if the CpG site is modified
    //         if (is_cpg)
    //         {
    //             // Increment the modified count
    //             modified_count++;
    //         } else {
    //             // Increment the unmodified count
    //             unmodified_count++;
    //         }
    //     }
    // }

    std::cout << "=== CpG Modification Counts ===" << std::endl;
    std::cout << "Modified CpG sites: " << modified_count << std::endl;
    std::cout << "Unmodified CpG sites: " << unmodified_count << std::endl;
    std::cout << "Total CpG sites: " << modified_count + unmodified_count << std::endl;
    std::cout << "Percentage of CpG sites modified: " << (double)modified_count / (modified_count + unmodified_count) * 100 << "%" << std::endl;
    return std::make_pair(modified_count, unmodified_count);
}

// Function to check if a given position is a CpG site (input is 1-indexed)
bool RefQuery::isCpG(std::string chr, int64_t pos, int strand)
{
    pos--;  // Convert to 0-indexed

    // Choose the map based on the strand
    std::map<std::string, std::map<int64_t, bool>>* map;
    if (strand == 0)
    {
        map = &this->chr_pos_to_cpg;
    } else {
        map = &this->chr_pos_to_cpg;
    }

    // Check if the chromosome is in the map
    if (map->find(chr) == map->end())
    {
        return false;
    }

    // Check if the position is in the map
    if ((*map)[chr].find(pos) == (*map)[chr].end())
    {
        return false;
    }

    // Return if the position is a CpG site
    return (*map)[chr][pos];

    // // Get the sequence
    // const std::string& sequence = this->chr_to_seq[chr];

    // // If on the forward strand, check if the position is a CpG site
    // if (strand == 0)
    // {
    //     if (pos + 1 < (int64_t)sequence.size() && sequence[pos] == 'C' && sequence[pos + 1] == 'G')
    //     {
    //         return true;
    //     }

    // // If on the reverse strand, check if the position is a CpG site
    // } else {
    //     if (pos - 1 >= 0 && sequence[pos] == 'G' && sequence[pos - 1] == 'C')
    //     {
    //         return true;
    //     }
    // }

    // return false;
    // Check if the position is in the map

    // 

    // // Determine the map to use based on the strand
    // std::map<std::string, std::map<int64_t, bool>>* map;
    // if (strand == 0)
    // {
    //     map = &this->chr_pos_to_cpg;
    // } else {
    //     map = &this->chr_pos_to_cpg_rev;
    // }

    // // Check if the chromosome is in the map
    // if (map->find(chr) == map->end())
    // {
    //     return false;
    // }

    // // Check if the position is in the map
    // if ((*map)[chr].find(pos) == (*map)[chr].end())
    // {
    //     return false;
    // }

    // // Return if the position is a CpG site
    // return (*map)[chr][pos];

    // // Check if the chromosome is in the map
    // if (this->chr_pos_to_cpg.find(chr) == this->chr_pos_to_cpg.end())
    // {
    //     return false;
    // }

    // // Check if the position is in the map
    // if (this->chr_pos_to_cpg[chr].find(pos) == this->chr_pos_to_cpg[chr].end())
    // {
    //     return false;
    // }

    // // Return if the position is a CpG site
    // return this->chr_pos_to_cpg[chr][pos];
}

int32_t RefQuery::getCpGSiteCount()
{
    return this->cpg_site_count;
}

std::string RefQuery::getFilepath()
{
    return this->fasta_filepath;
}

// Function to get the reference sequence at a given position range
char RefQuery::getBase(std::string chr, int64_t pos)
{
    // Convert positions from 1-indexed (reference) to 0-indexed (string indexing)
    pos--;

    // Ensure that the position is not negative
    if (pos < 0)
    {
        return 'N';
    }

    // Get the sequence
    const std::string& sequence = this->chr_to_seq[chr];

    // Get the base
    char base = sequence[pos];

    // If the base is empty, return empty string
    if (base == '\0')
    {
        return 'N';
    }

    return base;
}
