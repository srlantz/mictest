#include "HitStructures.h"

void BunchOfHits::SortByPhiBuildPhiBins()
{
  std::sort(&m_hits[0], &m_hits[m_fill_index], sortHitsByPhiMT);

  int last_bin = -1;
  int idx      =  0;
  for (int i = 0; i < m_fill_index; ++i)
  {
    Hit &h = m_hits[i];

    int bin = getPhiPartition(h.phi());

    if (bin != last_bin) 
    {
      m_phi_bin_infos[bin].first  = idx;
      // PhiBinInfo.second set to 0 in Reset()
    }
    ++m_phi_bin_infos[bin].second;

    last_bin = bin;
    ++idx;
  }

  //now fix empty bins
  int nextHitToFind = 0;
  for (int b = 0; b < Config::nPhiPart; ++b)
  {
    if (m_phi_bin_infos[b].first == 0 && m_phi_bin_infos[b].second == 0) 
    {
      m_phi_bin_infos[b].first = nextHitToFind;
    } 
    else 
    {
      nextHitToFind = m_phi_bin_infos[b].first + m_phi_bin_infos[b].second;
    }
    //std::cout << "bin=" << b << " set to " << m_phi_bin_infos[b].first << "," << m_phi_bin_infos[b].second << std::endl;
  }
}
