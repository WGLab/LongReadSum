#include "ComStruct.h"

Fast5ReaderRunOption::Fast5ReaderRunOption(){
   /*
   f5events.resize(MAX_F5EVENT_SIZE);
   f5eventsOld0.resize(MAX_F5EVENT_SIZE);
   f5eventsOld1.resize(MAX_F5EVENT_SIZE);
   f5signals.resize(MAX_F5SIGNAL_SIZE);
   f5moves.resize(MAX_F5SIGNAL_SIZE);
   group_dif.resize(MAX_F5SIGNAL_SIZE+1);
   group_sum.resize(MAX_F5SIGNAL_SIZE+1);
   */

   f5events = new F5Event[MAX_F5EVENT_SIZE];
   f5eventsOld0 = new F5EventOld0[MAX_F5EVENT_SIZE];
   f5eventsOld1 = new F5EventOld1[MAX_F5EVENT_SIZE];
   f5signals = new int16_t[MAX_F5SIGNAL_SIZE];
   f5moves = new int8_t[MAX_F5SIGNAL_SIZE];
   group_dif = new double[MAX_F5SIGNAL_SIZE+1];
   group_sum = new double[MAX_F5SIGNAL_SIZE+1];
}

Fast5ReaderRunOption::~Fast5ReaderRunOption(){
   delete [] f5events;
   delete [] f5eventsOld0;
   delete [] f5eventsOld1;
   delete [] f5signals;
   delete [] f5moves;
   delete [] group_dif;
   delete [] group_sum;
}


kmer_signal_model_struct_region::kmer_signal_model_struct_region(uint64_t p_ref_start_pos){
   region_kmer_signal.reserve(50000);
   ref_start_pos = p_ref_start_pos;
}
uint64_t kmer_signal_model_struct_region::get_ref_end_pos(){
   return ref_start_pos + region_kmer_signal.size();
}

uint64_t kmer_signal_model_struct_region::get_ref_start_pos(){
   return ref_start_pos;
}
      
void kmer_signal_model_struct_region::set_ref_start_pos(uint64_t p_ref_start_pos){
   ref_start_pos = p_ref_start_pos;
   region_kmer_signal.clear();
}
      
const kmer_signal_model_struct& kmer_signal_model_struct_region::operator[](size_t idx) const{
   return region_kmer_signal[idx];
}
      
void kmer_signal_model_struct_region::add(const kmer_signal_model_struct& ksms){
   region_kmer_signal.push_back(ksms);
}

uint64_t kmer_signal_model_struct_region::size(){
   return region_kmer_signal.size();
}

MapPos1Adj::MapPos1Adj(const Map1BasePos& mps, uint16_t p_adj){
   map_pos = mps;
   adjust = p_adj;
   _t_pos_dif = 100;
   _t_minus_dif = 100;
   _t_plus_dif = 100;

}
void MapPos1Adj::add_adjust(uint16_t p_adj){
   adjust |= p_adj;
}

int MapPos1Adj::get_adjust_num(){
   int adj_num = 0;
   if ((adjust&Pos_Match)>0){ adj_num+=1; }
   if ((adjust&Minus_Match)>0){ adj_num+=1; }
   if ((adjust&Plus_Match)>0){ adj_num+=1; }
   return adj_num;
}

ComparedPositionWithSignal::ComparedPositionWithSignal(){ 
   signal1.reserve(500); 
   signal2.reserve(500);
   signal1_seg_num.reserve(500);
   signal2_seg_num.reserve(500);
   reset();
}
void ComparedPositionWithSignal::reset(){
    signal1.clear();
    signal2.clear();
    signal1_seg_num.clear();
    signal2_seg_num.clear();

    signal_st.mean1 = 0;
    signal_st.mean2 = 0;
    signal_st.std1 = 0;
    signal_st.std2 = 0;

    signal_st.depth1 = 0;
    signal_st.depth2 = 0;

    signal_st.z_st = 0;
    signal_st.z_st_nb = 0;

    ref_pos.map_pos = 0;
    ref_pos.ref_base = 0;
}



