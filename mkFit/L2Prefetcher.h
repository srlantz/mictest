#ifndef L2Prefetcher_h
#define L2Prefetcher_h

#include "SideThread.h"

#include "HitStructures.h"


class L2Prefetcher : public SideThread<const BunchOfHits*>
{
public:
  L2Prefetcher(int cpuid=-1, int cpuid_st=-1)
  {
    SpawnSideThread(cpuid, cpuid_st);
  }

  ~L2Prefetcher()
  {
    // printf("L2Prefetcher::~L2Prefetcher will try to join the side thread now ...\n");

    JoinSideThread();
  }

  // virtual
  void DoWorkInSideThread(const BunchOfHits* boh)
  {
    // double t_0 = dtime();

    const char *pos = (const char *) & boh->m_hits[0];
    const char *end = (const char *) & boh->m_hits[boh->m_fill_index + 1];

    while (pos < end)
    {
      _mm_prefetch(pos, _MM_HINT_T1);
      pos += 64;
    }

    // double t_1 = dtime();
    // printf("Prefetched %d kB in %f s, %f %f\n",
    //        (end - (const char *) & boh->m_hits[0]) / 1024 + 1,
    //        t_1 - t_0, t_0, t_1);
  }
};

#endif
