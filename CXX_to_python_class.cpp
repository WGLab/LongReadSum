#include "CXX_to_python_class.h"

//// function for Output_FA
//
//
Output_FA::Output_FA(){
   read_length_count = new int64_t[MAX_READ_LENGTH];
   for(int _i_=0; _i_<MAX_READ_LENGTH; _i_++){
      read_length_count[ _i_ ] = 0;
   }
}

Output_FA::~Output_FA(){
   delete [] read_length_count;
}


//// function for Output_BAM
//
//
Output_BAM::Output_BAM(){
   for(int _i_=0; _i_<MAX_MAP_QUALITY; _i_++){
      map_quality_distribution[ _i_ ] = 0;
   }
   
   mapped_read_length_count = new int64_t[MAX_READ_LENGTH];
   for(int _i_=0; _i_<MAX_READ_LENGTH; _i_++){
      mapped_read_length_count[ _i_ ] = 0;
   }
}

Output_BAM::~Output_BAM(){
   delete [] mapped_read_length_count;
}



//// function for Output_F5
//
//
//
Output_F5::Output_F5(){
   passed_reads_list = new int64_t[MAX_READ_LENGTH];
   failed_reads_list = new int64_t[MAX_READ_LENGTH];
   for(int _i_=0; _i_<MAX_READ_LENGTH; _i_++){
      passed_reads_list[ _i_ ] = 0;
      failed_reads_list[ _i_ ] = 0;
   }
   for(int _i_=0; _i_<MAX_SIGNAL_VALUE; _i_++){
      signal_range[ _i_ ] = 0;
   }
}

Output_F5::~Output_F5(){
   delete [] passed_reads_list;
   delete [] failed_reads_list;
}



