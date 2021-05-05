#include "CXX_to_python_class.h"

Output_Info::Output_Info(){
   error_flag = 0;
   error_str = "";
}


//// function for Basic_Seq_Statistics
//
//
Basic_Seq_Statistics::Basic_Seq_Statistics(){
   read_length_count = new int64_t[MAX_READ_LENGTH];
   for(int _i_=0; _i_<MAX_READ_LENGTH; _i_++){
      read_length_count[ _i_ ] = ZeroDefault;
   }
   nx_read_length[9] = ZeroDefault;
}

Basic_Seq_Statistics::~Basic_Seq_Statistics(){
   delete [] read_length_count;
}


//// function for Basic_Seq_Quality_Statistics
//
//
Basic_Seq_Quality_Statistics::Basic_Seq_Quality_Statistics(){
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

Basic_Seq_Quality_Statistics::~Basic_Seq_Quality_Statistics(){
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
   
}

Output_BAM::~Output_BAM(){
}



//// function for Output_F5
//
//
//
Basic_F5_Statistics::Basic_F5_Statistics(){
   read_length_list = new int64_t[MAX_READ_LENGTH];
   for(int _i_=0; _i_<MAX_READ_LENGTH; _i_++){
      read_length_list[ _i_ ] = ZeroDefault;
   }
   for(int _i_=0; _i_<MAX_SIGNAL_VALUE; _i_++){
      signal_range[ _i_ ] = ZeroDefault;
   }
}

Basic_F5_Statistics::~Basic_F5_Statistics(){
   delete [] read_length_list;
}



