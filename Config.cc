// Include of Config.h not yet required

namespace Config
{
  int nTracks = 20000;
  int nEvents = 10;

  // Multi threading and Clone engine configuration
  int   numThreadsFinder = 1;
  
  // GPU computations
  int   numThreadsEvents = 1;
  int   numThreadsReorg = 1;

#ifdef __MIC__
  int   numThreadsSimulation = 60;
#else
  int   numThreadsSimulation = 12;
#endif

  bool  clonerUseSingleThread  = false;
  int   finderReportBestOutOfN = 1;

  int   numSeedsPerTask = 32;

  bool  useCMSGeom = false;
  bool  readCmsswSeeds = false;

  bool  cf_seeding  = false;
  bool  cf_fitting  = false;

  bool  super_debug = false;
}
