#include "ComFunction.h"

#include "glob.h"

#include <fstream>
#include <math.h>
#include <algorithm>

#include <sys/stat.h>
#include <unistd.h>

#include <iostream>


#include <iostream>

char m_cwd[FILENAME_MAX];
char mer5[5];

int get_cwd(){
  if (!GetCWD(m_cwd, sizeof(m_cwd))){
     return 1;
  }
  return 0;
}

bool isExist(const char * mpath){
  bool isexist = false;
  std::ifstream fin(mpath);
  if (fin && fin.good()){ isexist = true; }
  fin.close();
  return isexist;
}

/*#include<stdio.h>
bool check_exist_files(std::string & fn){
   FILE *file_open;
   if (file_open=open(fn.c_str(), "r")) {
       fclose(file_open);
       return true;
   }else{
       return false;
   }
}*/
bool check_exist_files(std::string & fn){
   if (access(fn.c_str(), F_OK) == -1){
      return false;
   }else{
      return true;
   }
}

int check_exist_folder(std::string & fn, bool create_model){
   struct stat folder_stats;
   int _stat_rt;
   //std::cout<<"check_exist_folder: for path="<<fn<<std::endl;
   _stat_rt = stat(fn.c_str(), &folder_stats);
   int is_exist=0;
   if (_stat_rt==0 && S_ISDIR(folder_stats.st_mode)){
      is_exist = 1;
      //std::cout<<"check_exist_folder_yes: for path="<<fn<<std::endl;
   }else{
      if (create_model){
          char m_c_path[256];
          char *move_p = NULL;
          size_t p_len;
          snprintf(m_c_path, sizeof(m_c_path),"%s", fn.c_str());
          p_len = fn.size(); //strlen(m_c_path);
          //std::cout<<"\tTo Create folder: "<< m_c_path<<" "<<p_len<<std::endl;
          if(m_c_path[p_len - 1] == '/') {  m_c_path[p_len - 1] = 0; }
          for(move_p = m_c_path + 1; *move_p; move_p++){
              if(*move_p == '/') {
                  *move_p = 0;
                  //std::cout<<"\tfor path="<<m_c_path<<std::endl;
                  struct stat folder_stats_in;
                  _stat_rt = stat(m_c_path, &folder_stats_in);
                  if (_stat_rt==0 && S_ISDIR(folder_stats_in.st_mode)){
                     //std::cout<<"\tfor path_existing="<<m_c_path<<" "<< _stat_rt <<std::endl; 
                     *move_p = '/';
                     continue;
                  }
                  //std::cout<<"\tfor path_Nexisting="<<m_c_path<< " " << _stat_rt<<std::endl;
                  //std::cout<<"create path"<<m_c_path<<std::endl;
                  is_exist = mkdir(m_c_path,0777);
                  //sleep(30);
                  if (is_exist==-1){
                     //std::cout<<"\tError Cannot create "<<fn<<std::endl;
                     *move_p = '/';
                     break;
                  }
                  //std::cout<<"\tCreate path "<<m_c_path<< " " <<is_exist <<" done"<<std::endl;
                  *move_p = '/';
                  //sleep(10);
              }
          }
          //std::cout<<"\tTo Final_path="<< m_c_path<<std::endl;
          struct stat folder_stats_fal;
          _stat_rt = stat(m_c_path, &folder_stats_fal);
          if (_stat_rt==0 && S_ISDIR(folder_stats_fal.st_mode)){ 
             //std::cout<<"\tTo Final_path_existing="<< m_c_path<<" "<< S_ISDIR(folder_stats_fal.st_mode) <<" "<< folder_stats_fal.st_mode <<" " << _stat_rt <<std::endl;
             is_exist = 0; 
          }else{
             //std::cout<<"create path"<<m_c_path<<std::endl;
             is_exist = mkdir(m_c_path,0777);
             //std::cout<<"\t create final_path"<<m_c_path<< " " <<is_exist << " " << _stat_rt <<" done"<<std::endl;
             //sleep(30);
          }
          /*if (mkdir(fn.c_str(), 0777)==-1){
             //std::cout<<"\tCannot create "<<fn<<std::endl;
             is_exist = -1;
          }else{
             is_exist = 0;
          }*/
      }else{
         is_exist = 0;
      }
   }
   return is_exist;
}



bool m_cp_str(char * dest, const char * msource, int mlen){
   int scp_res = snprintf(dest, mlen, "%s", msource);
   if (scp_res<1 || scp_res>=mlen){
       fprintf(stderr, "Error !!! Input (%s) is larger than the maximum (%d)", msource, mlen);
       return false;
   }
   return true;
}



int get_Files_From_Pattern(const std::string & file_pat, std::vector<std::string> * pat_files){
   glob_t glob_results;
   glob(file_pat.c_str(), GLOB_TILDE, NULL, &glob_results);
   for(unsigned int ig = 0; ig<glob_results.gl_pathc; ++ig){
      pat_files->push_back(std::string(glob_results.gl_pathv[ig]));
   }
   globfree(&glob_results);

   return 0;
}

int get_Files_From_Pattern(const std::string & file_pat, std::map<std::string, std::string> * file_dict ){
   std::map<std::string, std::string>::iterator fd_it;
   std::size_t slash_pos;

   glob_t glob_results;
   glob(file_pat.c_str(), GLOB_TILDE, NULL, &glob_results);
   for(unsigned int ig = 0; ig<glob_results.gl_pathc; ++ig){
      std::string cur_full_filename(glob_results.gl_pathv[ig]);
      slash_pos = cur_full_filename.find_last_of("/");
      std::string cur_fn;
      if (slash_pos==std::string::npos){
         cur_fn = cur_full_filename;
      }else{
         cur_fn = cur_full_filename.substr(slash_pos+1);
      }
      fd_it = file_dict->find(cur_fn);
      if (fd_it==file_dict->end()){
         (*file_dict)[cur_fn] = cur_full_filename;
      }else{
         fprintf(stdout, "Warning!!! duplicate files: %s \t %s\n", cur_full_filename.c_str(), (*file_dict)[cur_fn].c_str());
      }
   }
   globfree(&glob_results);

   return 0;
}

double st_mean(const uint64_t start_pos, const uint64_t length, const std::vector<double>& dvlist, const std::vector<int> ind_list){
    if (length==0){ return 0; }
    double v_sum = 0;
    int avail_num = 0;
    for (uint64_t dvi=0; dvi<length; dvi++){
        if (ind_list[start_pos+dvi]<1){ continue; }
        avail_num += 1;
        v_sum += dvlist[start_pos+dvi];
    }
    return (avail_num>0?v_sum/avail_num:0); 
} 

double st_std(const uint64_t start_pos, const uint64_t length, const double mean, const std::vector<double>& dvlist, const std::vector<int> ind_list){
    if (length==0){ return 0; }
    double std_sum = 0;
    int avail_num = 0;
    for (uint64_t dvi=0; dvi<length; dvi++){
        if (ind_list[start_pos+dvi]<1){ continue; }
        avail_num += 1;
        std_sum += pow(dvlist[start_pos+dvi]-mean, 2);
    }
    return (avail_num>1?sqrt(std_sum/(avail_num-1)):0);
}

double st_mean(const uint64_t start_pos, const uint64_t length, const std::vector<double>& dvlist){
   if (length==0){ return 0; }
   double v_sum = 0;
   for (uint64_t dvi=0; dvi<length; dvi++){
      //std::cout << dvlist[start_pos+dvi] << " ";
      v_sum += dvlist[start_pos+dvi];
   }
   //std::cout<< std::endl;
   return (length>0?v_sum/length:0);
}
double st_std(const uint64_t start_pos, const uint64_t length, const double mean, const std::vector<double>& dvlist){
   if (length==0){ return 0; }
   double std_sum = 0;
   //std::cout << " in std mean=" << mean << std::endl;
   for (uint64_t dvi=0; dvi<length; dvi++){
      std_sum += pow(dvlist[start_pos+dvi]-mean, 2);
      //std::cout << dvlist[start_pos+dvi] << " " << dvlist[start_pos+dvi]-mean << " " << pow(dvlist[start_pos+dvi]-mean, 2) << " " << std_sum << std::endl;
   }
   //return (length>0?sqrt(std_sum/(length)):0);
   return (length>1?sqrt(std_sum/(length-1)):0);
}

bool compare_RankPos_v (const RankPos& rp1, const RankPos& rp2) {
    if (rp1.value==rp2.value){ return (rp1.pos<rp2.pos); }
    else{ return (rp1.value<rp2.value); }
}

bool compare_RankPos_p (const RankPos& rp1, const RankPos& rp2) {
    return (rp1.pos<rp2.pos);
}

std::vector<RankPos> get_top_N_extreme(const double * data, const uint64_t m_start, uint64_t m_end, uint16_t topN, uint16_t sep_dist, int min_max){
   std::vector<RankPos> m_rank; m_rank.reserve(500);
   std::vector<RankPos> top_rank; top_rank.reserve(topN+5);

   if (m_end < sep_dist * 2){
      std::cout << "Error!!!! The end position is too small!!!" << m_end << "<" << sep_dist<<"*2" << std::endl;
      RankPos cur_rp;
      cur_rp.value = -1;
      cur_rp.pos = (m_start + m_end)/2;
      top_rank.push_back(cur_rp);
      return top_rank;
   }

   for(uint64_t pos_ind=(m_start>=sep_dist?m_start-sep_dist:0); pos_ind<m_end-sep_dist; pos_ind++){
       RankPos cur_rp;
       cur_rp.value = min_max>0?(data[pos_ind]):(-(data[pos_ind]));
       cur_rp.pos = pos_ind+sep_dist;
       m_rank.push_back(cur_rp);
   }
   
   std::sort(m_rank.begin(), m_rank.end(), compare_RankPos_v);
   /*
   for (int ssi=0; ssi<(m_rank.size()<50?m_rank.size():50); ssi++){
      std::cout << std::setprecision (5) << m_rank[ssi].value << "/" << m_rank[ssi].pos << " ";
   }
   std::cout<<std::endl; */

   uint64_t cal_sep_dist = (m_end - m_start)/(topN+1);
   if (cal_sep_dist > sep_dist ) {cal_sep_dist = sep_dist;}
   if (cal_sep_dist<3){
       /*std::cout<< "Warning!!! too less points for split: " << m_start << "-" << m_end << "=" << m_end-m_start << " sep_dist=" << sep_dist << "/(cal)" << cal_sep_dist<<" topN="<<topN;
       std::sort(m_rank.begin(), m_rank.end(), compare_RankPos_p);
       std::cout << ": ";
       for (size_t ssi=0; ssi<m_rank.size(); ssi++){
         std::cout << std::setprecision (3) << m_rank[ssi].value << "/" << (m_rank[ssi].pos)%100 << " ";
       }
       std::cout<<std::endl;*/
       for (int _c_i_t=0; _c_i_t<topN; _c_i_t++){
          if (m_start + (_c_i_t+1)*(cal_sep_dist>0?cal_sep_dist:1) < m_end - (cal_sep_dist>0?cal_sep_dist:1) ){
             RankPos cur_rp;
             cur_rp.value = -1;
             cur_rp.pos = m_start + (_c_i_t+1)*(cal_sep_dist>0?cal_sep_dist:1);
             top_rank.push_back(cur_rp);
          }
       }
   }else {
       std::vector<RankPos>::reverse_iterator mr_it;
       std::vector<RankPos>::iterator tr_it;
       for (mr_it=m_rank.rbegin(); mr_it!=m_rank.rend(); mr_it++){
          if (m_start>=sep_dist) {
             if (mr_it->pos - m_start < cal_sep_dist-1) { continue; }
          }
          if (m_end - mr_it->pos < cal_sep_dist) { continue; }

          for(tr_it=top_rank.begin(); tr_it!=top_rank.end(); tr_it++){
             if ( ((tr_it->pos>mr_it->pos)?(tr_it->pos - mr_it->pos):(mr_it->pos-tr_it->pos)) < cal_sep_dist){
                 break;
             }
          }
          if (tr_it==top_rank.end()){
             top_rank.push_back(*mr_it);
          }
          if (top_rank.size()>=topN){
             break;
          }
       }
   }
   /*if (!(cal_sep_dist<3) && top_rank.size() < topN){
      std::cout<< "Warning!!! less split points generated: " << m_start << "-" << m_end << "=" << m_end-m_start << " sep_dist=" << sep_dist << "/(cal)" << cal_sep_dist << " " << top_rank.size() << " " << topN;
       std::sort(m_rank.begin(), m_rank.end(), compare_RankPos_p);
       std::cout << ": ";
       for (size_t ssi=0; ssi<m_rank.size(); ssi++){
         for (std::vector<RankPos>::iterator c_rp_it=top_rank.begin(); c_rp_it!=top_rank.end(); c_rp_it++){
            if (m_rank[ssi].pos == c_rp_it->pos){
               std::cout<< " >< ";
               break;
            }
         }
         std::cout << std::setprecision (3) << m_rank[ssi].value << "/" << (m_rank[ssi].pos)%100 << " ";
       }
      std::cout<<std::endl;
   }*/
   std::sort(top_rank.begin(), top_rank.end(), compare_RankPos_p);

   return top_rank;
}

std::vector<std::string> m_split_string(const std::string & m_str, std::string multi_delimiters, bool contain_delimiter){
    std::vector<std::string> m_substr_list;
    m_substr_list.reserve(10);

    size_t cur_pos = 0;
    size_t last_pos = 0;
    while ((cur_pos=m_str.find_first_of(multi_delimiters, last_pos))!=std::string::npos){
        m_substr_list.push_back(m_str.substr(last_pos, cur_pos-last_pos+(contain_delimiter?1:0)));
        last_pos = cur_pos+1;
        if (last_pos==std::string::npos || last_pos>=m_str.size()){
            break;
        }
    }
    if (!(last_pos==std::string::npos || last_pos>m_str.size())){
        m_substr_list.push_back(m_str.substr(last_pos));
    }

    return m_substr_list;
}

std::string substr_to_lastN(const std::string p_str, int last2){
   std::string n_str3 = p_str;
   while (last2>0){
      std::size_t m_found = n_str3.find_last_of("/\\");
      if (std::string::npos == m_found){
          return n_str3;
      }
      n_str3 = n_str3.substr(0, m_found);
      last2--;
   }
   return n_str3;
}

std::string _remove_last_slash(std::string p_old_string){
   int slash_num = 0;
   while (p_old_string[p_old_string.size()-slash_num-1]=='/' || p_old_string[p_old_string.size()-slash_num-1]=='\\'){
      if (p_old_string.size()-slash_num > 1){
         slash_num++;
      }else{
          break;
      }
   }
   return p_old_string.substr(0, p_old_string.size()-slash_num);
}



