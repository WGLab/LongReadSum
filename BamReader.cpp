#include "BamReader.h"
#include <htslib/sam.h>

#include "ComFunction.h"

#include <iostream>
#include <sstream>
#include <fstream>

int BamReadOption::set_repdict1(const std::string & m_str, std::string m_delimiters){
    RepeatRegion crr;
    _c_oss_chr.str("");
    _c_oss_chr.clear();
    _c_oss_pos.str("");
    _c_oss_pos.clear();

    if (m_str.size()<3){
       return 0;
    }

    std::vector<std::string> rr_v = m_split_string(m_str, m_delimiters);
    //std::cout<<"\tTest split="<<rr_v[0]<<" "<<rr_v[1]<<" " <<rr_v[2]<<" "<<rr_v[3]<<" size="<<rr_v.size()<<std::endl;
    if (rr_v.size()<4){
       std::cout<< "\tWarning!!! Repeat pattern format is not correct! <" << m_str << ">" <<std::endl;
       return 1;
    }else{
       if (!m_cp_str(crr.chrn, rr_v[0].c_str(), CHAR_SIZE)){
          std::cout<< "\tWaring!!! Failed to get chr from " << m_str<< std::endl;
          return 3;
       }
       crr.start_pos = std::stoull(rr_v[1]);
       crr.end_pos = std::stoull(rr_v[2]);
       //crr.len_repeat_unit = std::stoi(rr_v[3]);
       crr.len_repeat_unit = rr_v[3].size();
    }
    //std::istringstream iss(m_str);
    //if (!(iss >> crr.chrn >> crr.start_pos >> crr.end_pos >> crr.len_repeat_unit)) { return 1; }
    if (crr.end_pos-crr.start_pos<1 || crr.len_repeat_unit<1){
        std::cout<<"\t Warning!!! Incorrect repeat info="<<m_str<<std::endl;
        return 2;
    }
    _c_oss_chr<<crr.chrn;
    _c_oss_pos<<crr.start_pos<<"-"<<crr.end_pos;
    //std::cout<< "read test " << crr.chrn << ":" << crr.start_pos<<"-"<<crr.end_pos << ":" << crr.len_repeat_unit<<std::endl;
    //std::cout<< "read test2 " << _c_oss_chr.str() << ":" << _c_oss_pos.str() <<std::endl;
    chr_dit = repdict.find(_c_oss_chr.str());
    if (chr_dit==repdict.end()){
       std::map<std::string, RepeatRegion> new_rr_dict;
       new_rr_dict[_c_oss_pos.str()] = crr;
       repdict[_c_oss_chr.str()] = new_rr_dict;
    }else{
       pos_it = chr_dit->second.find(_c_oss_pos.str());
       if (pos_it==chr_dit->second.end()){
          (chr_dit->second)[_c_oss_pos.str()] = crr;
       }else{
          std::cout<<"\t Warning!!! Repeat info already in the dict: "<<_c_oss_chr.str()<<":"<<_c_oss_pos.str()<< " ??? "<<chr_dit->first << ":" << pos_it->first<<std::endl;
      }
   }
   
   return 0;
}

int BamReadOption::set_repdict(std::string mrepfile){
   if (mrepfile.size()>0){
      //std::map<std::string, std::map<std::string, RepeatRegion> >::iterator chr_dit;
      //std::map<std::string, RepeatRegion>::iterator pos_it;
      //std::string line;
      //std::ostringstream oss_chr;
      //std::ostringstream oss_pos;

      //std::cout<<"Test repfile="<<mrepfile<<std::endl;
      std::ifstream infile(mrepfile);
      while (std::getline(infile, _c_line)){
         //std::cout<<"\tTest repline="<<_c_line<<std::endl;
         set_repdict1(_c_line);
          /*RepeatRegion crr;
          _c_oss_chr.clear();
          _c_oss_pos.clear();

          std::istringstream iss(_c_line);
          if (!(iss >> crr.chrn >> crr.start_pos >> crr.end_pos >> crr.len_repeat_unit)) { break; }
          if (crr.end_pos-crr.start_pos<1 || crr.len_repeat_unit<1){
              std::cout<<"\t Warning!!! Incorrect repeat info="<<line<<std::endl;
              continue;
          }
          _c_oss_chr<<crr.chrn;
          _c_oss_pos<<crr.start_pos<<"-"<<crr.end_pos;
          chr_dit = repdict.find(_c_oss_chr.str());
          if (chr_dit==repdict.end()){
             std::map<std::string, RepeatRegion> new_rr_dict;
             new_rr_dict[_c_oss_pos.str()] = crr;
             repdict[_c_oss_chr.str()] = new_rr_dict;
          }else{
             pos_it = chr_dit->second.find(_c_oss_pos.str());
             if (pos_it==chr_dit->second.end()){
                 (chr_dit->second)[_c_oss_pos.str()] = crr;
             }else{
                 std::cout<<"\t Warning!!! Repeat info already in the dict: "<<oss_chr.str()<<":"<<oss_pos.str()<< " ??? "<<chr_dit->first << ":" << pos_it->first<<std::endl;
             }
          }*/
      }
      infile.close();
   }
   return 0;
}

std::string BamReadOption::toString(){
   std::ostringstream oss_tostr;

   oss_tostr <<"\t" << (m_w_qry_seq?"With query sequence":"Without query sequence")<<"\n";
   oss_tostr <<"\t" << (m_w_qry_qual?"With query sequence quality":"Without query sequence quality")<<"\n";
   //oss_tostr <<"\t" << (m_w_map_detail?"With map detail":"Without map detail")<<"\n";
   oss_tostr <<"\t" << (m_w_pos_map_detail?"With map_pos detail":"Without map_pos detail")<<"\n";
   /*if (m_w_ref_seq){
      oss_tostr <<"\t" << "With ref sequence with the length=" << m_w_ref_seq_len<<"\n";
   }else{
      oss_tostr <<"\t" << "Without ref sequence"<<"\n";
   }*/
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

//int BamReadOption::set_w_map_detail(bool m_p_det, bool m_c_det, char * ref_seq, int m_refseq_len){
int BamReadOption::set_w_map_detail(bool m_p_det){
   //m_w_map_detail = m_c_det;
   m_w_pos_map_detail = m_p_det;
   //m_w_ref_seq = ref_seq;
   //m_w_ref_seq_len = m_refseq_len;
   return 0;
}

int BamReadOption::set_w_qry_qual(bool mqryqual){ m_w_qry_qual=mqryqual; return 0; }
bool BamReadOption::get_w_qry_qual(){ return m_w_qry_qual; }

bool BamReadOption::get_w_qry_seq(){ return m_w_qry_seq; }

//bool BamReadOption::get_w_map_detail(){return m_w_map_detail; }

bool BamReadOption::get_w_pos_map_detail(){ return m_w_pos_map_detail; }
/*char BamReadOption::get_ref_char_at_pos(uint64_t m_pos){
   if (m_w_map_detail) {
      std::cout<< "Warning!!! map_detail is not set" << std::endl;
      return 'N';
   }else if (m_pos>=m_w_ref_seq_len){
      std::cout<< "Warning!!! Pos=" << m_pos << " is larger than the length=" << m_w_ref_seq_len<< std::endl;
      return 'N';
   }else{
      return m_w_ref_seq[m_pos];
   }
}*/

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

void BamReadOption::check_specifiedRegion_multi(const char * chrn, uint64_t refStartPos, uint64_t refEndPos, std::vector<RepeatRegion>& c_rep_regions, std::map<std::string, std::map<std::string, RepeatRegion> > & p_rep_dict){
   chr_dit = p_rep_dict.find(chrn);
   if (chr_dit!=p_rep_dict.end()){
      for (pos_it=chr_dit->second.begin(); pos_it!=chr_dit->second.end(); pos_it++){
          ovlp_max_start = pos_it->second.start_pos > refStartPos ? pos_it->second.start_pos : refStartPos;
          ovlp_min_end = pos_it->second.end_pos < refEndPos ? pos_it->second.end_pos : refEndPos;
          if (ovlp_min_end >= ovlp_max_start + get_c_min_ovlp_len(refEndPos, refStartPos, pos_it->second.end_pos, pos_it->second.start_pos, min_ovlp_len)){
             RepeatRegion c_rr;
             if (m_cp_str(c_rr.chrn, chrn, CHAR_SIZE)){
                c_rr.start_pos = pos_it->second.start_pos;
                c_rr.end_pos = pos_it->second.end_pos;
                c_rr.len_repeat_unit = pos_it->second.len_repeat_unit;
                c_rep_regions.push_back(c_rr);
             } 
          }
      } 
   }
}
int BamReadOption::check_specifiedRegion(const char * chrn, uint64_t refStartPos, uint64_t refEndPos){
   if (m_w_specifiedRegion){
      //if (repdict.size()==0){
      //   return 2;
      //}

      chr_dit = repdict.find(chrn);
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


//////////////////////////////////////////////
const char BamReader::m_cigar_str[] = "MIDNSHP=XB";

int BamReader::openBam(const char * bamfile){
   //std::cout<<"In openBam1="<< bamfile<<" -> "<<_bam_file_str<<" "<< in_bam_file<<std::endl;
   if (!m_cp_str(in_bam_file, (char*)bamfile, CHAR_SIZE)) { return BAM_FAILED;}

   _bam_file_str = std::string(bamfile);

   bam_one_alignment = bam_init1();
   bam_1alignment_core = &bam_one_alignment->core;
  
   in_bam = sam_open(bamfile, "rb");
   if (NULL == in_bam){
      std::cout<< "Error! Cannot open sam file ("<<bamfile<< ")." << std::endl;
      bam_status = BAM_FAILED;
      init();
      return BAM_FAILED;
   }
   hdr = sam_hdr_read(in_bam);

   //std::cout<<"In openBam2="<< bamfile<<" -> "<<_bam_file_str<<" "<< in_bam_file<<std::endl;

   bam_status = BAM_OPEN;
   return 0;
}

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

int BamReader::resetBam(const char * bamfile){
   //std::cout<<"In resetBam="<< bamfile<<" -> "<<_bam_file_str<<" "<< in_bam_file<<std::endl;
   destroy();
   init();
   Bam1Record br1;
   br = br1;
   //br_list.clear();
   return openBam(bamfile);
}

//int BamReader::resetBam(){
//   return resetBam(_bam_file_str.c_str());
//   //return resetBam(in_bam_file);
//}

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

int BamReader::read1RecordFromBam(){
   if (!check_bam_status()){
      std::cout<< "No bam opened or Open bam failed." << std::endl;
      return 1;
   } 
 
   int failed = 2;
   //int sam_ret;

   //Bam1Record br1;
   //while (sam_ret = sam_read1(in_bam, hdr, bam_one_alignment)>=0){
   while (sam_read1(in_bam, hdr, bam_one_alignment)>=0){
      //if (_read1RecordFromBam_(br1)==0){
      if (_read1RecordFromBam_()==0){
          /*std::cout << "read1RecordFromBam" << "\n";
          std::cout << "Qry inf o= " << br.qry_name << ":" << br.qry_start_pos << "(" << br.qry_start_pos_rel << ")-" << br.qry_end_pos << " with qry seq length=" << br.qry_seq.size() <<"/" << br.qry_seq_len <<"\n";
          std::cout << "Map inf o= " << (br.map_strand==0?"+":"-") << br.map_chr << ":" << br.ref_start_pos << "-" << br.ref_end_pos <<"\n";
          std::cout << "Map flag =" << br.map_flag << " quality=" << br.map_qual << " cigar len=" << br.cigar_len.size() << "/" << br.cigar_type.size() <<"\n";*/
          failed = 0;
          break;
      }/*else{
          std::cout<<" Cannot read "<<Bam1Record_toString(br)<<std::endl;
       }*/
   }
 
   return failed;
}

int BamReader::reset_Bam1Record(){
   return reset_Bam1Record(br);
}
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
       //if (m_op == BAM_CSOFT_CLIP || BAM_CHARD_CLIP == m_op){
       //   std::cout<<" Debug " << br.qry_name << " " << mbp.qry_pos << "/" << m_len_original_read << " " <<mbp.ref_pos << " "<< m_op<<"/"<<mlen << std::endl;
       //}
   }
}

/*void BamReader::_set_map_detail(uint64_t ref_go_pos, uint64_t qry_go_pos, int mlen, Bam1Record & br, uint8_t m_op, bool ref_add, bool qry_add, const uint16_t m_map_strand, const uint64_t m_len_read){
   for (int m_movi=0; m_movi<mlen; m_movi++){
       Map1Base mb;
       mb.qry_base = (qry_add ? br.qry_seq[m_map_strand==0?(qry_go_pos+m_movi):(m_len_read - qry_go_pos - m_movi - 1)] : '-' );
       mb.ref_base = (ref_add ? bamReadOp.get_ref_char_at_pos(ref_go_pos + m_movi) : '-' );
       br.map_detail.push_back(mb);
   }
}
*/
int BamReader::_read1RecordFromBam_(){
   reset_Bam1Record(br);

   br.map_flag = bam_1alignment_core->flag;
   if (!bamReadOp.check_unmap(br.map_flag)){
      return 1;
   }
   if (br.map_flag & BAM_FUNMAP) {return 0;}

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
   br.map_chr = hdr->target_name[bam_1alignment_core->tid];
   
   br.qry_seq_len = bam_1alignment_core->l_qseq;
   if (!bamReadOp.check_min_qry_len(br.qry_seq_len)){
      return 4;
   }
   br.qry_name = bam_get_qname(bam_one_alignment);

   //br.qry_seq.clear();
   if (bamReadOp.get_w_qry_seq() || bamReadOp.get_w_pos_map_detail()){
     uint8_t * seq_int = bam_get_seq(bam_one_alignment);
     for (uint64_t  sqii=0; sqii < br.qry_seq_len; sqii++){
         br.qry_seq.push_back(seq_nt16_str[bam_seqi(seq_int, sqii)]);
      }
      // 012345678901234
      // ACMGRSVTWYHKDBN
      //
      // 1 for A, 2 for C, 4 for G,
      // 8 for T and 15 for N.
   }

   //br.qry_qual.clear(); 
   if (bamReadOp.get_w_qry_qual()){
      br.qry_qual.append(bam_get_qual(bam_one_alignment));
   }
 
   //br.map_qual = bam_1alignment_core->qual;
   //if ( br.map_qual < bamReadOp.get_map_quality_thr()){
   //   return 505;
   //}

   if(bamReadOp.get_w_specifiedRegion() && bamReadOp.check_specifiedRegion(br.map_chr.c_str(), br.ref_start_pos, br.ref_end_pos) < 0) {
      return 5;
   }

   // get qry information;
   uint32_t* m_cigar = bam_get_cigar(bam_one_alignment);
   //if (bamReadOp.get_w_unmap() && (BAM_FUNMAP & br.map_flag)){
   //   return 0;
   //}
   if (bamReadOp.get_w_pos_map_detail()){; // || bamReadOp.get_w_map_detail()){;
   }else{ return 0; }

   int m_op, m_next_op, m_pre_op;
   int m_len;

   br.left_clip = 0;
   br.right_clip = 0;
   ref_go_pos = br.ref_start_pos;
   qry_go_pos = 0; qry_go_pos_rel = 0;
   bool first_non_indel_clip = false;
   //uint64_t _len_original_read = 0;
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
             //std::cout<<" in read "<< ref_go_pos << " "<<  qry_go_pos << " " << m_len << std::endl;
             if (bamReadOp.get_w_pos_map_detail()) {_set_map_pos_detail(ref_go_pos, qry_go_pos, m_len, br, m_op, true, true, br.map_strand, br._len_original_read); }
             //if (bamReadOp.get_w_map_detail()) {_set_map_detail(ref_go_pos, qry_go_pos_rel, m_len, br, m_op, true, true, br.map_strand, _len_read); }
             //std::cout<<" in read "<< ref_go_pos << " "<<  qry_go_pos << " " << m_len << " " << br.map_pos_detail[br.map_pos_detail.size()-1].qry_pos << "/" << br.map_pos_detail[br.map_pos_detail.size()-1].ref_pos << std::endl;
             //
             //revised on Dec 17, 2020
             //br.qry_end_pos = qry_go_pos;
             ref_go_pos += m_len;
             qry_go_pos += m_len;
             qry_go_pos_rel += m_len;
             //revised on Dec 17, 2020
             br.qry_end_pos = qry_go_pos;
             break;
          case BAM_CINS:  // I  //fprintf(stdout, "%s%d", "I", m_len);
             if (bamReadOp.get_w_pos_map_detail()) { _set_map_pos_detail(ref_go_pos, qry_go_pos, m_len, br, m_op, false, true, br.map_strand, br._len_original_read); }
             //if (bamReadOp.get_w_map_detail()) { _set_map_detail(ref_go_pos, qry_go_pos_rel, m_len, br, m_op, false, true, br.map_strand, _len_read); }
             qry_go_pos += m_len;
             qry_go_pos_rel += m_len;
             break;
          case BAM_CDEL:   // D //fprintf(stdout, "%s%d", "D", m_len);
             if (bamReadOp.get_w_pos_map_detail()){ _set_map_pos_detail(ref_go_pos, qry_go_pos, m_len, br, m_op, true, false, br.map_strand, br._len_original_read); }
             //if (bamReadOp.get_w_map_detail()){ _set_map_detail(ref_go_pos, qry_go_pos_rel, m_len, br, m_op, true, false, br.map_strand, _len_read); }
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
             //if (bamReadOp.get_w_map_detail()) { _set_map_detail(ref_go_pos, qry_go_pos_rel, m_len, br, m_op, false, true, br.map_strand, _len_read); }
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
             //if (bamReadOp.get_w_map_detail()) { _set_map_detail(ref_go_pos, qry_go_pos_rel, m_len, br, m_op, true, true, br.map_strand, _len_read); }
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
    
   /*std::cout << "_read1RecordFromBam_" << "\n";
   std::cout << "Qry inf o= " << br.qry_name << ":" << br.qry_start_pos << "(" << br.qry_start_pos_rel << ")-" << br.qry_end_pos << " with qry seq length=" << br.qry_seq.size() <<"/" << br.qry_seq_len <<"\n";
   std::cout << "Map inf o= " << (br.map_strand==0?"+":"-") << br.map_chr << ":" << br.ref_start_pos << "-" << br.ref_end_pos <<"\n";
   std::cout << "Map flag =" << br.map_flag << " quality=" << br.map_qual << " cigar len=" << br.cigar_len.size() << "/" << br.cigar_type.size() <<"\n";
   */
   return 0;
}

int BamReader::readBam(std::vector<Bam1Record> &br_list, int num_records, bool br_clear){
   //std::cout<<bamReadOp.toString();

   if (!check_bam_status()){
      std::cout<< "No bam opened or Open bam failed." << std::endl;
      return 1;
   }
   if ( br_clear){
      br_list.clear();
   }
   //int sam_ret;
   t_num_records = 0;
   //while (sam_ret = sam_read1(in_bam, hdr, bam_one_alignment)>=0){
   while (sam_read1(in_bam, hdr, bam_one_alignment)>=0){
       //Bam1Record br1;
       //if (_read1RecordFromBam_(br1)==0){
       if (_read1RecordFromBam_()==0){
          //Bam1Record brc = br;
          // br_list.push_back(brc);
          br_list.push_back(br);
          t_num_records += 1;
          if (num_records>0 && t_num_records>=num_records){
             break;
          }
       }/*else{
          std::cout<<" Cannot read "<<Bam1Record_toString(br)<<std::endl;
       }*/
   }

   return 0;
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
   //oss_tostr << " Map-flag=" << br.map_flag << " qual=" << br.map_qual << " cigar-len=" << br.cigar_len.size() << "/" << br.cigar_type.size();
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
   return Bam1Record_toString(br);
}




