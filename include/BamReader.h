#ifndef BAMREADER_H_
#define BAMREADER_H_

#include <htslib/sam.h>
#include <htslib/bgzf.h>

#include <iostream>
#include <sstream>
#include <fstream>

#include <mutex>

#include "ComStruct.h"
#include "OutputStructures.h"

#define BAM_UN_OPEN 1
#define BAM_FAILED 2
#define BAM_OPEN 3
#define BAM_CLOSE 4

#define DNASEQ 1
#define RNASEQ 2

// Class representing an HTSlib BAM file
class HTSReader {
    public:
        htsFile* bam_file; // open the BAM file for reading
        bam_hdr_t* header; // read the BAM header
        bam1_t* record;
        int record_count = 0;

        // Bool for whether the reading is complete
        bool reading_complete = false;

        // Read the next batch of records from the BAM file
        int readNextRecords(int batch_size, Output_BAM & output_data);

        // Return if the reader has finished reading the BAM file
        bool hasNextRecord();

        HTSReader(const std::string &bam_file_name);
        ~HTSReader();
};

// Type representing data for a single record in the BAM file
typedef struct Bam1Record{
    int num_mismatch = -1;
    uint64_t qry_start_pos = 0;
    uint64_t qry_end_pos = 0;
    uint64_t qry_start_pos_rel = 0;
    uint64_t qry_seq_len = 0;
    uint64_t _len_original_read = 0;
    uc8string qry_qual;
    std::string qry_name;
    std::string qry_seq;

    uint16_t pri_sec_sup;
    uint16_t map_flag;
    uint8_t  map_qual;
    uint16_t  map_strand;
    std::vector<uint32_t> cigar_len;
    std::vector<uint16_t> cigar_type;

    std::vector<Map1Base> map_detail;
    std::vector<Map1BasePos> map_pos_detail;

    uint64_t ref_start_pos;
    uint64_t ref_end_pos;
    uint64_t left_clip;
    uint64_t right_clip;
    std::string map_chr;
} Bam1Record;

class BamReadOption{

private:
   std::map<std::string, std::map<std::string, RepeatRegion> >::iterator chr_dit;
   std::map<std::string, RepeatRegion>::iterator pos_it;
   bool ovlp_region;
   uint64_t ovlp_max_start;
   uint64_t ovlp_min_end;

   std::ostringstream _c_oss_chr;
   std::ostringstream _c_oss_pos;
   std::string _c_line;

   bool m_w_qry_seq;
   bool m_w_qry_qual;

   //bool m_w_map_detail;
   bool m_w_pos_map_detail;
   //char * m_w_ref_seq;
   //uint64_t m_w_ref_seq_len;

   bool m_w_unmap;
   bool m_w_supplementary;
   bool m_w_secondary;

   bool m_w_specifiedRegion;
   uint64_t min_ovlp_len;

   uint64_t min_read_len;

   uint8_t map_quality_thr;

public:
   uint8_t data_type; 
   uint8_t get_map_quality_thr();
   int set_map_quality_thr(uint8_t p_map_quality_thr);
   int set_w_qry_seq(bool mqryseq);
   int set_w_map_detail(bool m_p_det); //, bool m_c_det, char * ref_seq, int m_refseq_len);
   int set_w_qry_qual(bool mqryqual);
   bool get_w_qry_qual(); 
   bool get_w_qry_seq();
   //bool get_w_map_detail();
   bool get_w_pos_map_detail();
   //char get_ref_char_at_pos(uint64_t m_pos);
   uint64_t get_min_read_len(){return min_read_len;}
   std::map<std::string, std::map<std::string, RepeatRegion> > repdict;

   int set_min_read_len(uint64_t _min_read_len){ min_read_len = _min_read_len; return 0; }
   int set_w_unmap(bool m_unmap);
   int set_w_supplementary(bool m_sup);
   int set_w_secondary(bool m_sec);
   int set_w_specifiedRegion(bool spR, uint64_t mol);
   int set_w_specifiedRegion(bool spR);

   int set_repdict(std::string mrepfile);
   int set_repdict1(const std::string & m_str, std::string m_delimiter=WhiteSpace);
   //std::map<std::string, std::map<std::string, RepeatRegion> > get_repdict(){ return repdict; }

   bool get_w_unmap();
   bool get_w_supplementary();
   bool get_w_secondary();
   bool get_w_specifiedRegion();
   uint64_t get_min_ovlp_len(); 

   bool check_unmap(uint8_t curflag);
   bool check_supplementary(uint8_t curflag);
   bool check_secondary(uint8_t curflag);
   int check_specifiedRegion(const char * chrn, uint64_t refStartPos, uint64_t refEdnPos);
   void check_specifiedRegion_multi(const char * chrn, uint64_t refStartPos, uint64_t refEndPos, std::vector<RepeatRegion>& c_rep_regions, std::map<std::string, std::map<std::string, RepeatRegion> > & p_rep_dict);
   // not safe
   std::map<std::string, std::map<std::string, RepeatRegion> >::iterator get_chr_it();
   std::map<std::string, RepeatRegion>::iterator get_pos_it();
   bool check_min_qry_len(uint64_t _qry_len);

   BamReadOption();
   ~BamReadOption();

   std::string toString();
};

class BamReader{
    private:
        void init();
        void destroy();

        int bam_status;

        char in_bam_file[CHAR_SIZE];
        samFile * in_bam;
        bam_hdr_t * hdr;
        bam1_t * bam_one_alignment;
        bam1_core_t * bam_1alignment_core;
        uint64_t ref_go_pos;
        uint64_t qry_go_pos;
        uint64_t qry_go_pos_rel;

        void _set_map_pos_detail(uint64_t ref_go_pos, uint64_t qry_go_pos, int mlen, Bam1Record & br, uint8_t m_op, bool ref_add, bool qry_add, const uint16_t m_map_strand, const uint64_t m_len_original_read);
        //void _set_map_detail(uint64_t ref_go_pos, uint64_t qry_go_pos, int mlen, Bam1Record & br, uint8_t m_op, bool ref_add, bool qry_add, const uint16_t m_map_strand, const uint64_t m_len_read);

        int t_num_records;

        std::string _bam_file_str;

    public:
        std::vector<Bam1Record> record_list;  // Vector of individual BAM file records
        std::vector<Bam1Record> readNextNRecords(int n, std::vector<Bam1Record> & record_list);
        int readNextRecord(Bam1Record& br);

        const char * m_cigar_str = "MIDNSHP=XB";
        BamReadOption bamReadOp;
        int openBam(const char * bamfile);
        bool check_bam_status();

        Bam1Record current_bam_record;

        // Get the number of records in the bam file
        int getNumberOfRecords(const char * bamfile);
        int resetBam(const char * bamfile);
        //int resetBam();
        std::string get_bam_file();
        int reset_Bam1Record();
        int reset_Bam1Record(Bam1Record & br);

        // Get the record list
        std::vector<Bam1Record> getRecordList();

        static std::string Bam1Record_toString(Bam1Record & br);
        std::string Bam1Record_toString();
        static std::string Basic_Bam1Record_toString(Bam1Record & br);
        static std::string Min_Bam1Record_toString(Bam1Record & br);
        BamReader(char * bamfile);
        BamReader();
        ~BamReader();
};

#endif


