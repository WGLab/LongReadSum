def count_cpg_sites(genome_file):
    cpg_count = 0

    # Loop through each sequence in the genome file
    with open(genome_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                continue
            genome = line.upper()
            for i in range(len(genome) - 1):
                if genome[i:i+2] == 'CG':
                    cpg_count += 1

    # with open(genome_file, 'r') as file:
    #     genome = file.read().upper()

    # for i in range(len(genome) - 1):
    #     if genome[i:i+2] == 'CG':
    #         cpg_count += 1

    return cpg_count

# Usage example
genome_file = "/mnt/isilon/wang_lab/perdomoj/data/ont_aws/gm24385_q20_2021.10/GRCh38_no_alt_analysis_set.fa"
cpg_count = count_cpg_sites(genome_file)
print(f"Number of CpG sites: {cpg_count}")
