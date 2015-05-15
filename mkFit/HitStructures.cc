#include "HitStructures.h"

void EventOfCandidates::FillTrackPacks(std::vector<std::pair<int, int>>& seed_packs)
{
  // XXXX This is half cooked
  // Rethink what it should actually do.
  // Should we make sure there is no seeds from different eta bins in
  // the same pack?
  // Should we sort seeds in phi / pt (or some other metric that
  // takes road extent into account?

  int count = 0;
  for (auto &i : m_etabins_of_candidates)
  {
    count += i.m_fill_index / NN + 1;
  }

  seed_packs.resize(count);
  for (auto &t : seed_packs)
  {
  }

  for (auto &i : m_etabins_of_candidates)
  {
    int n = i.m_candidates.m_fill_index;
    for (int i = 0;n > 0)
    {
      
    }
  }
}
