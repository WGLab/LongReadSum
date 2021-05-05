#ifndef COMFUNCTION_H_
#define COMFUNCTION_H_

#include <stdio.h>

#include <map>
#include <vector>
#include <string>

#include <limits>
#include <iomanip>

#include "ComStruct.h"

#define get_array_size(m_a) sizeof(m_a)/sizeof(m_a[0])

#define round1(d) (((double)((int)(d*10+0.5)))/10)
#define round2(d) (((double)((int)(d*100+0.5)))/100)
#define round3(d) (((double)((int)(d*1000+0.5)))/1000)

#define UNUSED(expr) do { (void)(expr); } while (0)

#ifdef WINDOWS
    #include <direct.h>
    #define GetCWD _getcwd
#else
    #include <unistd.h>
    #define GetCWD getcwd
#endif


#define rev_base(b) (b=='A'?'T':(b=='C'?'G':(b=='G'?'C':(b=='T'?'A':(b=='U'?'A':'N')))))

#define ind_dna_base(b) (b=='A'?1:(b=='C'?2:(b=='G'?3:(b=='T'?4:(0)))))
#define ind_rna_base(b) (b=='A'?1:(b=='C'?2:(b=='G'?3:(b=='U'?4:(0)))))


extern char mer5[];



bool isExist(const char *);

bool check_exist_files(std::string & fn);
int check_exist_folder(std::string & fn, bool create_model);

int get_Files_From_Pattern(const std::string & file_pat, std::map<std::string, std::string> * file_dict);
int get_Files_From_Pattern(const std::string & file_pat, std::vector<std::string> * pat_files);

extern char m_cwd[];
int get_cwd();

double st_mean(const uint64_t start_pos, const uint64_t length, const std::vector<double>& dvlist);
double st_std(const uint64_t start_pos, const uint64_t length, const double mean, const std::vector<double>& dvlist);

double st_mean(const uint64_t start_pos, const uint64_t length, const std::vector<double>& dvlist, const std::vector<int> ind_list);
double st_std(const uint64_t start_pos, const uint64_t length, const double mean, const std::vector<double>& dvlist, const std::vector<int> ind_list);

bool compare_RankPos_v (const RankPos& rp1, const RankPos& rp2);
bool compare_RankPos_p (const RankPos& rp1, const RankPos& rp2); 

std::vector<RankPos> get_top_N_extreme(const double * data, const uint64_t m_start, uint64_t m_end, uint16_t topN=1, uint16_t sep_dist=4, int min_max=1);

bool m_cp_str(char * dest, const char * msource, int mlen);

std::vector<std::string> m_split_string(const std::string & m_str, std::string multi_delimiters=WhiteSpace, bool contain_delimiter=false);

std::string _remove_last_slash(std::string p_old_string);
std::string substr_to_lastN(const std::string p_str, int last2 = 2);

#endif
