#ifndef COMSTRUCT_H_
#define COMSTRUCT_H_

#include <string>
#include <map>
#include <vector>


/////////////////////////////////
#define MAX_F5EVENT_SIZE   3000000
#define MAX_F5SIGNAL_SIZE 15000000

#define KMER_SIZE 6

#define NORM_SIGNAL_RANGE 4
#define INVAIL_SIGNAL -1000

//////////////////////////////////
#define MODULE_NUM 5
#define CHAR_SIZE 1024

#define WhiteSpace "\t\n\v\f\r "

// for feature generation for deepmod
//                            ,  ,
// training need large          ,  ,
//#define DP_Max_Instance_size 1000000
//for testing with smaller      ,
#define DP_Max_Instance_size 100000
#define RNN_Window 21
#define Feature_Size 7

//////////////////////////////////

typedef struct F5Event{
   float mean;
   float stdv;
   uint64_t start;
   uint64_t length;
   char model_state[KMER_SIZE];
   uint32_t move;
} F5Event;

typedef struct F5EventOld1{
   float mean;
   float stdv;
   float start;
   float length;
   char model_state[KMER_SIZE];
   uint32_t move;
} F5EventOld1;

typedef struct F5EventOld0{
   double mean;
   double start;
   double stdv;
   double length;
   char model_state[KMER_SIZE];
   uint64_t move;
} F5EventOld0;


typedef struct F5AnnoEvent{
   float mean;
   float stdv;
   uint64_t start;
   uint64_t length;
   char model_state[KMER_SIZE];
   uint64_t ref_pos;
   uint64_t qry_pos;
   uint64_t map_type;
} F5AnnoEvent;

///////////////////////////////////////
typedef struct RankPos{
   double value;
   uint64_t pos;
} RankPos;


///////////////////////////////////
typedef std::basic_string<unsigned char> uc8string;

typedef struct Map1Base{
  char qry_base;
  char ref_base;
} Map1Base;
// 0-based, like bam and bed format
typedef struct Map1BasePos{
  uint64_t qry_pos;
  uint64_t ref_pos;
  uint64_t map_type;
} Map1BasePos;

#define Pos_Match 1
#define Minus_Match 2
#define Plus_Match 4
class MapPos1Adj{
  public:
    uint16_t adjust;
    Map1BasePos map_pos;

    double _t_pos_dif;
    double _t_minus_dif;
    double _t_plus_dif;

    MapPos1Adj(const Map1BasePos& mps, uint16_t p_adj=0);
    void add_adjust(uint16_t p_adj);
    int get_adjust_num();
};

typedef struct Map1BasePosPred{
  uint64_t qry_pos;
  uint64_t ref_pos;
  uint64_t map_type;
  
  //float pred1;
  //float pred2;
  float pred;
} Map1BasePosPred;

typedef struct MapRecord{
   std::string ref_chr;
   uint16_t ref_strand;
   std::string qry_readname;
   uint64_t _start_;

   //std::vector<Map1BasePos> map_detail;
   std::vector<Map1BasePosPred> map_detail;
} MapRecord;


////////////////////////////

/*typedef struct Fast5ReaderRunOption{
  uint64_t group_size; 
  std::string read_num;
  std::string read_id;
  F5Event f5events[MAX_F5EVENT_SIZE];
  F5EventOld0 f5eventsOld0[MAX_F5EVENT_SIZE];
  F5EventOld1 f5eventsOld1[MAX_F5EVENT_SIZE];
  int16_t f5signals [MAX_F5SIGNAL_SIZE];
  int8_t f5moves [MAX_F5SIGNAL_SIZE];
  double group_dif[MAX_F5SIGNAL_SIZE+1];
  double group_sum[MAX_F5SIGNAL_SIZE+1];
  //std::map<std::string, std::vector<std::string> > errorInfo;

  //std::map<std::string, std::string> readToF5;
  //std::map<std::string, bool> readOfInterest;
} Fast5ReaderRunOption;
*/

class Fast5ReaderRunOption{
public:
  uint64_t group_size;
  std::string read_num;
  std::string read_id;
  /*std::vector<F5Event> f5events;
  std::vector<F5EventOld0> f5eventsOld0;
  std::vector<F5EventOld1> f5eventsOld1;
  std::vector<int16_t> f5signals;
  std::vector<int8_t> f5moves;
  std::vector<double> group_dif;
  std::vector<double> group_sum;*/

  F5Event *f5events;
  F5EventOld0 *f5eventsOld0;
  F5EventOld1 *f5eventsOld1;
  int16_t * f5signals;
  int8_t * f5moves;
  double *group_dif;
  double * group_sum;

  //F5Event f5events[MAX_F5EVENT_SIZE];
  //F5EventOld0 f5eventsOld0[MAX_F5EVENT_SIZE];
  //F5EventOld1 f5eventsOld1[MAX_F5EVENT_SIZE];
  //int16_t f5signals [MAX_F5SIGNAL_SIZE];
  //int8_t f5moves [MAX_F5SIGNAL_SIZE];
  //double group_dif[MAX_F5SIGNAL_SIZE+1];
  //double group_sum[MAX_F5SIGNAL_SIZE+1];

  Fast5ReaderRunOption();
  ~Fast5ReaderRunOption();
};

// 0-based, like bam and bed format; end-not-included 
typedef struct GenomicRegion{
   uint64_t start_pos;
   uint64_t end_pos;
   char chrn[CHAR_SIZE];
} GenomicRegion;

// 0-based, like bam and bed format; end-not-included
typedef struct RepeatRegion{
   uint64_t start_pos;
   uint64_t end_pos;
   char chrn[CHAR_SIZE];
   int len_repeat_unit;
} RepeatRegion;

//////////////////////////////

struct kmer_signal_model_struct{
   float signal_mean;
   float signal_std;
};

class kmer_signal_model_struct_region {
   private:
      std::vector<kmer_signal_model_struct> region_kmer_signal;
      uint64_t ref_start_pos;
   public:
      kmer_signal_model_struct_region(uint64_t p_ref_start_pos);
      uint64_t get_ref_end_pos();
      uint64_t get_ref_start_pos();
      void set_ref_start_pos(uint64_t p_ref_start_pos);
      const kmer_signal_model_struct& operator[](size_t idx) const;
      void add(const kmer_signal_model_struct& ksms);
      uint64_t size();
};

struct Ref_Position_Base{
   uint64_t map_pos;
   char ref_base;
};

struct Ref_Position{
   uint64_t _position_;
   uint16_t _strand_;
};

struct Signal_ST{
   double z_st;
   double mean1;
   double mean2;
   double std1;
   double std2;
  
   uint64_t depth1;
   uint64_t depth2;

   double z_st_nb;
};

class ComparedPositionWithSignal{
   public: 
      Signal_ST signal_st;
      Ref_Position_Base ref_pos;
      std::vector<double> signal1;
      std::vector<double> signal2;
      std::vector<int> signal1_seg_num;
      std::vector<int> signal2_seg_num;
      
      ComparedPositionWithSignal(); 
      void reset();
};

struct F5AnnoIndexRecord{
   std::string f5_file;
   std::string qry_name;
   std::string f5anno_file;
   uint16_t pri_sec_sup;

   //uint16_t map_strand;
   std::string map_strand;
   std::string map_chr;
   uint64_t ref_start_pos;
   uint64_t ref_end_pos;
};


struct Pred_Mod_Info{
   char strand;
   uint64_t ref_pos;
   
   uint64_t mod_coverage;
   uint64_t unmod_coverage;
};


#endif

