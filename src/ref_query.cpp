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
                    this->chr_pos_to_cpg[chr][pos] = true;
                }
            }
        }
    }
    std::cout << "CpG sites located." << std::endl;
}

// Function to check if a given position is a CpG site (input is 1-indexed)
bool RefQuery::isCpG(std::string chr, int64_t pos)
{
    pos--;  // Convert to 0-indexed
    
    // Check if the chromosome is in the map
    if (this->chr_pos_to_cpg.find(chr) == this->chr_pos_to_cpg.end())
    {
        return false;
    }

    // Check if the position is in the map
    if (this->chr_pos_to_cpg[chr].find(pos) == this->chr_pos_to_cpg[chr].end())
    {
        return false;
    }

    // Return if the position is a CpG site
    return this->chr_pos_to_cpg[chr][pos];
}
