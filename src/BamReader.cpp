/*

BamReader.cpp:
Class for reading a set number of records from a BAM file. Used for multi-threading.

*/


#include <iostream>
#include <sstream>
#include <fstream>
#include <htslib/sam.h>

#include "BamReader.h"
#include "ComFunction.h"

int BamReadOption::set_repdict1(const std::string & input_pattern, std::string str_delimiters){
    RepeatRegion repeat_region;  // Where the repeat region will be stored
    _c_oss_chr.str("");
    _c_oss_chr.clear();
    _c_oss_pos.str("");
    _c_oss_pos.clear();

    if (input_pattern.size()<3){
       return 0;
    }

    std::vector<std::string> region_splitstr = readRepeatDataFromString(input_pattern, str_delimiters);
    // TODO: The following tests should occur within the function
    if (region_splitstr.size()<4){
       std::cout<< "\tWarning!!! Repeat pattern format is not correct! <" << input_pattern << ">" <<std::endl;
       return 1;
    }else{
       if (!isExpectedLength(repeat_region.repeat_size, region_splitstr[0].c_str(), CHAR_SIZE)){
          std::cout<< "\tUnexpected repeat pattern length:  " << input_pattern<< std::endl;
          return 3;
       }
       repeat_region.start_pos = std::stoull(region_splitstr[1]);
       repeat_region.end_pos = std::stoull(region_splitstr[2]);
       repeat_region.len_repeat_unit = region_splitstr[3].size();
    }
    if (repeat_region.end_pos-repeat_region.start_pos<1 || repeat_region.len_repeat_unit<1){
        std::cout<<"\t Warning!!! Incorrect repeat info="<<input_pattern<<std::endl;
        return 2;
    }
    _c_oss_chr<<repeat_region.repeat_size;
    _c_oss_pos<<repeat_region.start_pos<<"-"<<repeat_region.end_pos;
    chr_dit = repdict.find(_c_oss_chr.str());
    if (chr_dit==repdict.end()){
       std::map<std::string, RepeatRegion> new_rr_dict;
       new_rr_dict[_c_oss_pos.str()] = repeat_region;
       repdict[_c_oss_chr.str()] = new_rr_dict;
    }else{
       pos_it = chr_dit->second.find(_c_oss_pos.str());
       if (pos_it==chr_dit->second.end()){
          (chr_dit->second)[_c_oss_pos.str()] = repeat_region;
       }else{
          std::cout<<"\t Warning!!! Repeat info already in the dict: "<<_c_oss_chr.str()<<":"<<_c_oss_pos.str()<< " ??? "<<chr_dit->first << ":" << pos_it->first<<std::endl;
      }
   }
   
   return 0;
}

int BamReadOption::set_repdict(std::string mrepfile){
   if (mrepfile.size()>0){
      std::ifstream infile(mrepfile);
      while (std::getline(infile, _c_line)){
         set_repdict1(_c_line);
      }
      infile.close();
   }
   return 0;
}

std::string BamReadOption::toString(){
   std::ostringstream oss_tostr;

   oss_tostr <<"\t" << (m_w_qry_seq?"With query sequence":"Without query sequence")<<"\n";
   oss_tostr <<"\t" << (m_w_qry_qual?"With query sequence quality":"Without query sequence quality")<<"\n";
   oss_tostr <<"\t" << (m_w_pos_map_detail?"With map_pos detail":"Without map_pos detail")<<"\n";
   oss_tostr <<"\t" << (m_w_unmap?"With unmapped reads":"Without unmapped reads")<<"\n";
   oss_tostr <<"\t" << (m_w_supplementary?"With supplementary map":"Without supplementary map")<<"\n";
   oss_tostr <<"\t" << (m_w_secondary?"With seconday map":"Without seconday map")<<"\n";
   if (m_w_specifiedRegion){
       oss_tostr <<"\t" << "Specify a region of interest and require "<<min_ovlp_len<< " bp overlap"<<"\n";
       for (chr_dit=repdict.begin(); chr_dit!=repdict.end(); chr_dit++){
           for (pos_it=chr_dit->second.begin(); pos_it!=chr_dit->second.end(); pos_it++){
               oss_tostr <<"\t\t" << (chr_dit->first) << ":"<<(pos_it->first) << ":" << pos_it->second.len_repeat_unit << "\n";
           }
       }
   }else{
       oss_tostr <<"\t" << "Get all alignment rather than for specific regions" <<"\n";
   }

   return oss_tostr.str();
}

int BamReadOption::set_w_qry_seq(bool mqryseq){ m_w_qry_seq = mqryseq; return 0; }

int BamReadOption::set_w_map_detail(bool m_p_det){
   m_w_pos_map_detail = m_p_det;
   return 0;
}

int BamReadOption::set_w_qry_qual(bool mqryqual){ m_w_qry_qual=mqryqual; return 0; }
bool BamReadOption::get_w_qry_qual(){ return m_w_qry_qual; }

bool BamReadOption::get_w_qry_seq(){ return m_w_qry_seq; }

bool BamReadOption::get_w_pos_map_detail(){ return m_w_pos_map_detail; }

BamReadOption::BamReadOption(){
    data_type = DNASEQ;
    set_w_qry_seq(false);
    set_w_map_detail(true); //, false, NULL, 0);
    set_w_unmap(false);
    set_w_supplementary(true);
    set_w_secondary(true);
    set_w_specifiedRegion(false);
    set_w_qry_qual(false);
    set_min_read_len(50);
    set_map_quality_thr(0);
}

BamReadOption::~BamReadOption(){ set_w_map_detail(true); }//, false, NULL, 0); }

int BamReadOption::set_w_unmap(bool m_unmap){ m_w_unmap = m_unmap; return 0; }
int BamReadOption::set_w_supplementary(bool m_sup){ m_w_supplementary = m_sup; return 0; }
int BamReadOption::set_w_secondary(bool m_sec){ m_w_secondary = m_sec; return 0; }
int BamReadOption::set_w_specifiedRegion(bool spR){ return set_w_specifiedRegion(spR, 10); }
int BamReadOption::set_w_specifiedRegion(bool spR, uint64_t mol){ 
   m_w_specifiedRegion = spR; 
   min_ovlp_len = mol;
   if (mol<1 && m_w_specifiedRegion){
      std::cout<<"Warning!!! the length for overlaping (" << mol << ") is too small."<< std::endl;
   }
   return 0;
}
uint8_t BamReadOption::get_map_quality_thr(){
   return map_quality_thr;
}
int BamReadOption::set_map_quality_thr(uint8_t p_map_quality_thr){
   map_quality_thr = p_map_quality_thr;
   return 0;
}

bool BamReadOption::get_w_unmap(){ return m_w_unmap; }
bool BamReadOption::get_w_supplementary(){ return m_w_supplementary; }
bool BamReadOption::get_w_secondary(){ return m_w_secondary; }
bool BamReadOption::get_w_specifiedRegion(){ return m_w_specifiedRegion; }
uint64_t BamReadOption::get_min_ovlp_len(){ return min_ovlp_len; }

bool BamReadOption::check_unmap(uint8_t curflag){
   if ((!m_w_unmap) && (curflag & BAM_FUNMAP)) { return false; }
   else {return true;}
}
bool BamReadOption::check_supplementary(uint8_t curflag){
   if ((!m_w_supplementary) && (curflag & BAM_FSUPPLEMENTARY)) { return false; }
   else {return true;}
}
bool BamReadOption::check_secondary(uint8_t curflag){
   if ((!m_w_secondary) && (curflag & BAM_FSECONDARY)) { return false; }
   else {return true;}
}

bool BamReadOption::check_min_qry_len(uint64_t _qry_len){
   if (_qry_len > min_read_len){ return true; }
   else { return false; }
}

inline uint64_t get_c_min_ovlp_len(uint64_t endp1, uint64_t startp1, uint64_t endp2, uint64_t startp2, uint64_t p_min_ovlp_len){
   uint64_t d1 = (endp1 > startp1+5) ? (endp1 - startp1) : 5;
   uint64_t d2 = (endp2 > startp2+5) ? (endp2 - startp2) : 5;
   uint64_t mmin = d1>d2 ? d2 : d1;
   mmin = mmin > p_min_ovlp_len ? p_min_ovlp_len : mmin;
   if (mmin < 5) { return 5; }
   else { return mmin; }
}

void BamReadOption::check_specifiedRegion_multi(const char * repeat_size, uint64_t refStartPos, uint64_t refEndPos, std::vector<RepeatRegion>& all_repeat_regions, std::map<std::string, std::map<std::string, RepeatRegion> > & p_rep_dict){
   chr_dit = p_rep_dict.find(repeat_size);
   if (chr_dit!=p_rep_dict.end()){
      for (pos_it=chr_dit->second.begin(); pos_it!=chr_dit->second.end(); pos_it++){
          ovlp_max_start = pos_it->second.start_pos > refStartPos ? pos_it->second.start_pos : refStartPos;
          ovlp_min_end = pos_it->second.end_pos < refEndPos ? pos_it->second.end_pos : refEndPos;
          if (ovlp_min_end >= ovlp_max_start + get_c_min_ovlp_len(refEndPos, refStartPos, pos_it->second.end_pos, pos_it->second.start_pos, min_ovlp_len)){
             RepeatRegion repeat_region;
             if (isExpectedLength(repeat_region.repeat_size, repeat_size, CHAR_SIZE)){
                repeat_region.start_pos = pos_it->second.start_pos;
                repeat_region.end_pos = pos_it->second.end_pos;
                repeat_region.len_repeat_unit = pos_it->second.len_repeat_unit;
                all_repeat_regions.push_back(repeat_region);
             } 
          }
      } 
   }
}
int BamReadOption::check_specifiedRegion(const char * repeat_size, uint64_t refStartPos, uint64_t refEndPos){
   if (m_w_specifiedRegion){
      chr_dit = repdict.find(repeat_size);
      if (chr_dit==repdict.end()){
         return -1;
      }

      ovlp_region = false;
      for (pos_it=chr_dit->second.begin(); pos_it!=chr_dit->second.end(); pos_it++){
          ovlp_max_start = pos_it->second.start_pos > refStartPos ? pos_it->second.start_pos : refStartPos;
          ovlp_min_end = pos_it->second.end_pos < refEndPos ? pos_it->second.end_pos : refEndPos;
          if (ovlp_min_end >= ovlp_max_start + get_c_min_ovlp_len(refEndPos, refStartPos, pos_it->second.end_pos, pos_it->second.start_pos, min_ovlp_len)){
              ovlp_region = true;
              break;
          }
      }
      if (ovlp_region){ return 2; }
      else { return -2;}
   }else{
      return 1;
   }
}

std::map<std::string, std::map<std::string, RepeatRegion> >::iterator BamReadOption::get_chr_it(){
   return chr_dit;
}
std::map<std::string, RepeatRegion>::iterator BamReadOption::get_pos_it(){
   return pos_it;
}


const char BamReader::m_cigar_str[] = "MIDNSHP=XB";

int BamReader::getNumberOfRecords(const char * bamfile){
    htsFile* bam_file = hts_open(bamfile, "r"); // open the BAM file for reading
    bam_hdr_t* header = sam_hdr_read(bam_file); // read the BAM header
    bam1_t* record = bam_init1(); // allocate memory for a single record

    int record_count = 0;
    while (sam_read1(bam_file, header, record) >= 0) {
        // process the record here, if needed
        record_count++;
    }

    bam_destroy1(record); // free the memory for the record
    bam_hdr_destroy(header); // free the memory for the header
    hts_close(bam_file); // close the BAM file

    std::cout << "Number of records: " << record_count << std::endl;
    return record_count;
}

int BamReader::openBam(const char * bamfile){
   if (!isExpectedLength(in_bam_file, (char*)bamfile, CHAR_SIZE)) { return BAM_FAILED;}

   _bam_file_str = std::string(bamfile);

   bam_one_alignment = bam_init1();
   bam_1alignment_core = &bam_one_alignment->core;
  
   in_bam = sam_open(bamfile, "rb");
   if (NULL == in_bam){
      std::cerr<< "Cannot open SAM file ("<<bamfile<< ")." << std::endl;
      bam_status = BAM_FAILED;
      init();
      return BAM_FAILED;
   }
   hdr = sam_hdr_read(in_bam);
   bam_status = BAM_OPEN;

   return 0;
}

// Create a new bam reader with the given bam file
BamReader::BamReader(const char * bamfile){
   openBam(bamfile);
}

void BamReader::init(){
   in_bam = NULL;
   hdr = NULL;
   bam_one_alignment = NULL;
   bam_1alignment_core = NULL;
   bam_status = BAM_UN_OPEN;
}

BamReader::BamReader(){
   in_bam_file[0] = '\0';
   init();
}

// Destroy the bam reader and create a new one with the given bam file
int BamReader::resetBam(const char * bamfile){
   destroy();
   init();
   Bam1Record br1;
   this->current_bam_record = br1;
   //br_list.clear();
   return openBam(bamfile);
}

std::string BamReader::get_bam_file(){
   return _bam_file_str;
}

BamReader::~BamReader(){
   destroy();
}
void BamReader::destroy(){
   if (bam_one_alignment!=NULL){
      bam_destroy1(bam_one_alignment);
   }
   if (hdr!=NULL){
      bam_hdr_destroy(hdr);
   }
   if (in_bam!=NULL){
      sam_close(in_bam);
   }
   bam_status = BAM_CLOSE;
   init();
}

bool BamReader::check_bam_status(){
   if (bam_status==BAM_OPEN){
      return true;
   }else{
      return false;
   }
}

// Read one record from the bam file
int BamReader::read1RecordFromBam(){
   if (!check_bam_status()){
      std::cerr << "Unable to read BAM record" << std::endl;
      return 1;
   } 
 
   int failed = 2;
   while (sam_read1(in_bam, hdr, bam_one_alignment)>=0){
      if (_read1RecordFromBam_()==0){
          failed = 0;
          break;
      }
   }
 
   return failed;
}

int BamReader::reset_Bam1Record(){
   return reset_Bam1Record(this->current_bam_record);
}

// Clear the BAM record structure
int BamReader::reset_Bam1Record(Bam1Record & br){
   br.cigar_len.clear();
   br.cigar_type.clear();
   if (br.cigar_len.capacity()<1000){ br.cigar_len.reserve(1000); }
   if (br.cigar_type.capacity()<1000){ br.cigar_type.reserve(1000); }

   br.map_detail.clear();
   br.map_pos_detail.clear();
   if (br.map_pos_detail.capacity()<50000){ br.map_pos_detail.reserve(50000); }

   br.qry_qual.clear();
   br.qry_name.clear();
   br.qry_seq.clear();

   br.map_chr.clear();

   br.qry_seq.clear();
   br.qry_qual.clear();

   br._len_original_read = 0;

   return 0;
}

void BamReader::_set_map_pos_detail(uint64_t ref_go_pos, uint64_t qry_go_pos, int mlen, Bam1Record & br, uint8_t m_op, bool ref_add, bool qry_add, const uint16_t m_map_strand, const uint64_t m_len_original_read){
   //std::cout<<" in _set_map_pos_detail "<<ref_go_pos<< " "<<  qry_go_pos << " " << mlen << std::endl;
   uint64_t opsize = mlen;
   opsize = opsize << 4;
   for (int m_movi=0; m_movi<mlen; m_movi++){
       Map1BasePos mbp;
       mbp.qry_pos = qry_go_pos + (qry_add?m_movi:0);
       if (m_map_strand!=0){
          if (mbp.qry_pos>=m_len_original_read){
              std::cout<< "Error!!! mbp.qry_pos("<< mbp.qry_pos << ") >= m_len_original_read(" <<m_len_original_read<<")"<<std::endl;
          }
          mbp.qry_pos = m_len_original_read - mbp.qry_pos - 1; 
       }
       mbp.ref_pos = ref_go_pos + (ref_add?m_movi:0);
       mbp.map_type = m_op;
       mbp.map_type = mbp.map_type | opsize;
       br.map_pos_detail.push_back(mbp);
   }
}

int BamReader::_read1RecordFromBam_(){
    // Update the BAM record with a new instance
    Bam1Record br;
    this->current_bam_record = br;

   // Get the number of mismatched bases using the MD tag
   uint8_t *nmTag = bam_aux_get(bam_one_alignment, "NM");
   if (nmTag != NULL) {
      br.num_mismatch = bam_aux2i(nmTag);
   } else {
      br.num_mismatch = -1;
   }

   br.qry_name = bam_get_qname(bam_one_alignment);
   br.map_flag = bam_1alignment_core->flag;
   if (!bamReadOp.check_unmap(br.map_flag)){
      return 1;
   }
   if (!bamReadOp.check_supplementary(br.map_flag)){
      return 2;
   }
   if (!bamReadOp.check_secondary(br.map_flag)){
      return 3;
   }
   br.map_qual = bam_1alignment_core->qual;
   if ( br.map_qual < bamReadOp.get_map_quality_thr()){
      return 505;
   }

   br.pri_sec_sup = 0;
   if (br.map_flag & BAM_FSUPPLEMENTARY) { br.pri_sec_sup = 1000; }
   if (br.map_flag & BAM_FSECONDARY) { br.pri_sec_sup = 100; }
   br.ref_start_pos = bam_1alignment_core->pos;
   br.ref_end_pos = bam_endpos(bam_one_alignment);

   br.map_strand = 0; // forward;
   if (bam_is_rev(bam_one_alignment) ) { br.map_strand = 1; }
   if ( !( br.map_flag & BAM_FUNMAP ) ){
      br.map_chr = hdr->target_name[bam_1alignment_core->tid];
   }else{ br.map_chr = ""; }

   br.qry_seq_len = bam_1alignment_core->l_qseq;  // Get the sequence length
   if ( br.map_flag & BAM_FSECONDARY ){
      size_t cal_len = bam_cigar2qlen(bam_1alignment_core->n_cigar, bam_get_cigar(bam_one_alignment));
      if ( cal_len > br.qry_seq_len){ br.qry_seq_len = cal_len; }
   }
   if (!bamReadOp.check_min_qry_len(br.qry_seq_len)){
      return 4;
   }
   if ( (bamReadOp.get_w_qry_seq() || bamReadOp.get_w_pos_map_detail()) && (!( br.map_flag & BAM_FSECONDARY )) ){
     uint8_t * seq_int = bam_get_seq(bam_one_alignment);
     for (uint64_t  sqii=0; sqii < br.qry_seq_len; sqii++){
         br.qry_seq.push_back(seq_nt16_str[bam_seqi(seq_int, sqii)]);
      }
   }

   if (bamReadOp.get_w_qry_qual()){
      if (!( br.map_flag & BAM_FSECONDARY )) { br.qry_qual.append(bam_get_qual(bam_one_alignment)); }
   }

   if(bamReadOp.get_w_specifiedRegion() && bamReadOp.check_specifiedRegion(br.map_chr.c_str(), br.ref_start_pos, br.ref_end_pos) < 0) {
      return 5;
   }
   // get qry information;
   uint32_t* m_cigar = bam_get_cigar(bam_one_alignment);
   
   int m_op, m_next_op, m_pre_op;
   int m_len;
   if (bamReadOp.get_w_pos_map_detail()){; // || bamReadOp.get_w_map_detail()){;
   }else{ 
      for (uint32_t  m_i_cigar=0; m_i_cigar<bam_1alignment_core->n_cigar; ++m_i_cigar){
         m_op = bam_cigar_op(m_cigar[m_i_cigar]); br.cigar_type.push_back(m_op);
         m_len = bam_cigar_oplen(m_cigar[m_i_cigar]); br.cigar_len.push_back(m_len);
      }
      return 0; 
   }

   br.left_clip = 0;
   br.right_clip = 0;
   ref_go_pos = br.ref_start_pos;
   qry_go_pos = 0; qry_go_pos_rel = 0;
   bool first_non_indel_clip = false;
   br._len_original_read = 0;
   uint64_t _len_read = 0;
   uint64_t _len_align = 0;
   for (uint32_t m_i_cigar=0; m_i_cigar<bam_1alignment_core->n_cigar; ++m_i_cigar){
       m_op = bam_cigar_op(m_cigar[m_i_cigar]); //br.cigar_type.push_back(m_op);
       if (m_op==BAM_CDEL || m_op==BAM_CREF_SKIP || m_op==BAM_CPAD || m_op==BAM_CBACK){
          continue;
       }
       br._len_original_read += bam_cigar_oplen(m_cigar[m_i_cigar]); //br.cigar_len.push_back(m_len);
       if (m_op!=BAM_CHARD_CLIP){
          _len_read += bam_cigar_oplen(m_cigar[m_i_cigar]); //br.cigar_len.push_back(m_len);
       }
       if ((m_op!=BAM_CHARD_CLIP) && (m_op!=BAM_CSOFT_CLIP)){
          _len_align += bam_cigar_oplen(m_cigar[m_i_cigar]);
       }       
   }

   if (bamReadOp.get_min_read_len()>_len_align){return 6;}

   for (uint32_t  m_i_cigar=0; m_i_cigar<bam_1alignment_core->n_cigar; ++m_i_cigar){ 
      m_op = bam_cigar_op(m_cigar[m_i_cigar]); br.cigar_type.push_back(m_op);
      m_len = bam_cigar_oplen(m_cigar[m_i_cigar]); br.cigar_len.push_back(m_len);

      if (m_i_cigar==0){
         if (m_op==BAM_CSOFT_CLIP || m_op==BAM_CHARD_CLIP){
            if (br.map_strand==0) { br.left_clip = m_len; }
            else{ br.right_clip = m_len;}
            if (m_i_cigar+1<bam_1alignment_core->n_cigar){
               m_next_op = bam_cigar_op(m_cigar[m_i_cigar+1]);
               if (m_next_op!=BAM_CEQUAL && m_next_op!=BAM_CMATCH){
                  std::cout<<"Warning-non-match! The second aligned bases are not matached "<< br.qry_name <<" in " << _bam_file_str <<" " << m_i_cigar+1<<":"<< m_next_op << "/" << bam_cigar_oplen(m_cigar[m_i_cigar+1]) <<std::endl;
               }
            }
         }else if (m_op!=BAM_CEQUAL && m_op!=BAM_CMATCH){
            std::cout<<"Warning-non-match! The first aligned bases are not matached "<< br.qry_name<<" in " << _bam_file_str <<" " << m_i_cigar<<":"<< m_op <<"/"<< m_len <<std::endl;
         }
      }else if (m_i_cigar == bam_1alignment_core->n_cigar-1){
         if (m_op==BAM_CSOFT_CLIP || m_op==BAM_CHARD_CLIP){
             if (br.map_strand==0) { br.right_clip = m_len; }
             else{ br.left_clip = m_len; }
             if (m_i_cigar>0){
                 m_pre_op = bam_cigar_op(m_cigar[m_i_cigar-1]);
                 if (m_pre_op !=BAM_CEQUAL && m_pre_op!=BAM_CMATCH){
                     std::cout<<"Warning-non-match! The last second aligned bases are not matached "<< br.qry_name<<" in " << _bam_file_str <<" " << m_i_cigar-1<<":"<< m_pre_op << "/" << bam_cigar_oplen(m_cigar[m_i_cigar-1]) <<std::endl;
                 }
             }
         }else if (m_op!=BAM_CEQUAL && m_op!=BAM_CMATCH){
            std::cout<<"Warning-non-match! The last aligned bases are not matached "<< br.qry_name <<" in " << _bam_file_str <<" " << m_i_cigar<<":"<< m_op<<"/"<< m_len <<std::endl;
         }
      }

      switch (m_op){
          case BAM_CEQUAL:
          case BAM_CMATCH: // M
             if (!first_non_indel_clip){
                 br.qry_start_pos = qry_go_pos;
                 br.qry_start_pos_rel = qry_go_pos_rel;
                 first_non_indel_clip = true;
             }
             if (bamReadOp.get_w_pos_map_detail()) {_set_map_pos_detail(ref_go_pos, qry_go_pos, m_len, br, m_op, true, true, br.map_strand, br._len_original_read); }
             ref_go_pos += m_len;
             qry_go_pos += m_len;
             qry_go_pos_rel += m_len;
             br.qry_end_pos = qry_go_pos;
             break;
          case BAM_CINS: // I
             if (bamReadOp.get_w_pos_map_detail()) { _set_map_pos_detail(ref_go_pos, qry_go_pos, m_len, br, m_op, false, true, br.map_strand, br._len_original_read); }
             qry_go_pos += m_len;
             qry_go_pos_rel += m_len;
             break;
          case BAM_CDEL:   // D //fprintf(stdout, "%s%d", "D", m_len);
             if (bamReadOp.get_w_pos_map_detail()){ _set_map_pos_detail(ref_go_pos, qry_go_pos, m_len, br, m_op, true, false, br.map_strand, br._len_original_read); }
             ref_go_pos += m_len;
             break;
          case BAM_CREF_SKIP: // R //fprintf(stdout, "%s%d", "R", m_len);
             //_set_map_pos_detail(ref_go_pos, qry_go_pos, m_len, br, m_op, true, false);
             if (bamReadOp.data_type==DNASEQ){
                fprintf(stdout, "Caution for DNA data: N (skip) cigar %d:%d exists.\n", m_op, m_len);
             }
             ref_go_pos += m_len;
             break;
          case BAM_CSOFT_CLIP:  // S //fprintf(stdout, "%s%d", "S", m_len);
             if (bamReadOp.get_w_pos_map_detail()){ _set_map_pos_detail(ref_go_pos, qry_go_pos, m_len, br, m_op, false, true, br.map_strand, br._len_original_read); }
             qry_go_pos += m_len;
             qry_go_pos_rel += m_len;
             break; 
          case BAM_CHARD_CLIP: // H //fprintf(stdout, "%s%d", "H", m_len);
             //for getting end position in reads;
             if (bamReadOp.get_w_pos_map_detail()){ _set_map_pos_detail(ref_go_pos, qry_go_pos, m_len, br, m_op, false, true, br.map_strand, br._len_original_read); }
             ////if (bamReadOp.get_w_map_detail()) { _set_map_detail(ref_go_pos, qry_go_pos_rel, m_len, br, m_op, false, true, br.map_strand, _len_read); }
             qry_go_pos += m_len;
             break;
          case BAM_CPAD:  // P //fprintf(stdout, "%s%d", "P", m_len);
             break;
          case BAM_CDIFF: // X /fprintf(stdout, "%s%d", "X", m_len);
             if (bamReadOp.get_w_pos_map_detail()) { _set_map_pos_detail(ref_go_pos, qry_go_pos, m_len, br, m_op, true, true, br.map_strand, br._len_original_read); }
             ref_go_pos += m_len;
             qry_go_pos += m_len;
             qry_go_pos_rel += m_len;
             break;
          default:
             fprintf(stderr, "Unknow cigar %d:%d\n", m_op, m_len);
      }   
   }
  
   if (br.qry_end_pos>0){ // to include br.qry_end_pos
      br.qry_end_pos = br.qry_end_pos - 1;
   }else{
      return -1;
   }
   return 0;
}

// Read a batch of records from the bam file
//int BamReader::readBam(std::vector<Bam1Record> &record_list, int batch_size, bool clear_records){
int BamReader::readBam(){
    if (!check_bam_status()){
        std::cout<< "No bam opened or Open bam failed." << std::endl;
        return 1;
    }

   // Clear the record list if a new file is opened
    this->record_list.clear();
//   if (clear_records){
////      record_list.clear();
//   }

    this->t_num_records = 0;  // Reset the number of records
    while (sam_read1(in_bam, hdr, bam_one_alignment)>=0){
        // sam_read1 returns 0 if success
        if ((_read1RecordFromBam_())==0){
            this->record_list.push_back(this->current_bam_record);
            t_num_records += 1;

//            // If the number of records reaches the batch size, break
//            if (batch_size > 0 && t_num_records >= batch_size){
//                break;
//            }
        }
   }

   return 0;
}

// Get the list of records
std::vector<Bam1Record> BamReader::getRecordList() {
   return this->record_list;
}

std::string BamReader::Min_Bam1Record_toString(Bam1Record & br){
  std::ostringstream oss_tostr;
  oss_tostr << br.qry_name << ":" << br.qry_start_pos<<"-" << br.qry_end_pos<<"/"<<br._len_original_read;
  oss_tostr << " to " << (br.map_strand==0?"+":"-") << br.map_chr << ":" << br.ref_start_pos << "-" << br.ref_end_pos ;
  return oss_tostr.str();
}

std::string BamReader::Basic_Bam1Record_toString(Bam1Record & br){
   std::ostringstream oss_tostr;

   oss_tostr << "Qry= " << br.qry_name << ":" << br.qry_start_pos << "(" << br.qry_start_pos_rel << ")-" << br.qry_end_pos << " with=" << br.qry_seq.size() <<"/" << br.qry_seq_len<<"/"<<br._len_original_read ;
   
   oss_tostr << " Map= " << (br.map_strand==0?"+":"-") << br.map_chr << ":" << br.ref_start_pos << "-" << br.ref_end_pos ;
   oss_tostr << " Map-flag=" << br.map_flag <<                            " cigar-len=" << br.cigar_len.size() << "/" << br.cigar_type.size();

   return oss_tostr.str();
}


std::string BamReader::Bam1Record_toString(Bam1Record & br){
   std::ostringstream oss_tostr; 

   size_t topn = 3, topi=0;

   oss_tostr << "Qry inf o= " << br.qry_name << ":" << br.qry_start_pos << "(" << br.qry_start_pos_rel << ")-" << br.qry_end_pos << " with qry seq length=" << br.qry_seq.size() <<"/" << br.qry_seq_len <<"\n";
   oss_tostr << "Map inf o= " << (br.map_strand==0?"+":"-") << br.map_chr << ":" << br.ref_start_pos << "-" << br.ref_end_pos <<"\n";
   oss_tostr << "Map flag =" << br.map_flag << " quality=" << +(br.map_qual) << " cigar-len=" << br.cigar_len.size() << "/" << br.cigar_type.size() <<"\n";
   oss_tostr << "   ";
   for (size_t cgi=0; cgi<(br.cigar_len.size()>topn?topn:br.cigar_len.size()); cgi++){
      oss_tostr << br.cigar_type[cgi]<<":"<<br.cigar_len[cgi] << " ";
   }
   oss_tostr <<"\n";
   std::vector<Map1BasePos>::iterator mbp1_it;
   oss_tostr <<"Map_pos=" << br.map_pos_detail.size() << "/"<< (br.map_pos_detail.size() - br.left_clip - br.right_clip) << "\n";
   for (mbp1_it=br.map_pos_detail.begin(); mbp1_it!=br.map_pos_detail.end(); mbp1_it++){
      if (topi>=topn){break;}
      oss_tostr << "\t" <<mbp1_it->ref_pos<<"<-->"<<mbp1_it->qry_pos<<" <---"<<mbp1_it->map_type;
      topi++;
   }
   oss_tostr << "\n";

   return oss_tostr.str();
}

std::string BamReader::Bam1Record_toString(){
   return Bam1Record_toString(this->current_bam_record);
}




