#include "CXX_to_python_class.h"

//// function for Output_FA
//
//
Output_FA::Output_FA(){
   read_length_count = new int64_t[MAX_READ_LENGTH];
   for(int _i_=0; _i_<MAX_READ_LENGTH; _i_++){
      read_length_count[ _i_ ] = ZeroDefault;
   }
   nx_read_length[9] = ZeroDefault;
}

Output_FA::~Output_FA(){
   delete [] read_length_count;
}


//// function for Output_BAM
//
//
Output_FQ::Output_FQ(){
   pos_quality_distribution = new int[MAX_READ_LENGTH];
   pos_quality_distribution_dev = new float[MAX_READ_LENGTH];
   pos_quality_distribution_count = new int64_t[MAX_READ_LENGTH];
   for(int _i_=0; _i_<MAX_READ_LENGTH; _i_++){
      pos_quality_distribution[ _i_ ] = ZeroDefault;
      pos_quality_distribution_dev[ _i_ ] = ZeroDefault;
      pos_quality_distribution_count[ _i_ ] = ZeroDefault;
   }
   for(int _i_=0; _i_<MAX_BASE_QUALITY; _i_++){
      base_quality_distribution[ _i_ ] = ZeroDefault;
   }
}

Output_FQ::~Output_FQ(){
   delete [] pos_quality_distribution;
   delete [] pos_quality_distribution_dev;
   delete [] pos_quality_distribution_count;
}

//// function for Output_BAM
//
//
Output_BAM::Output_BAM(){
   for(int _i_=0; _i_<MAX_MAP_QUALITY; _i_++){
      map_quality_distribution[ _i_ ] = ZeroDefault;
   }
   
   mapped_read_length_count = new int64_t[MAX_READ_LENGTH];
   for(int _i_=0; _i_<MAX_READ_LENGTH; _i_++){
      mapped_read_length_count[ _i_ ] = ZeroDefault;
   }
   mapped_nx_read_length[9] = ZeroDefault;
}

Output_BAM::~Output_BAM(){
   delete [] mapped_read_length_count;
}



//// function for Output_F5
//
//
//
Output_F5::Output_F5(){
   passed_read_length_list = new int64_t[MAX_READ_LENGTH];
   failed_read_length_list = new int64_t[MAX_READ_LENGTH];
   for(int _i_=0; _i_<MAX_READ_LENGTH; _i_++){
      passed_read_length_list[ _i_ ] = ZeroDefault;
      failed_read_length_list[ _i_ ] = ZeroDefault;
   }
   for(int _i_=0; _i_<MAX_SIGNAL_VALUE; _i_++){
      signal_range[ _i_ ] = ZeroDefault;
   }
}

Output_F5::~Output_F5(){
   delete [] passed_read_length_list;
   delete [] failed_read_length_list;
}



