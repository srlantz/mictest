#include "Track.h"

#if 0
void Track::write_out(FILE *fp)
{
  Track t = clone_for_io();
  fwrite(&t, sizeof(Track), 1, fp);

  int nh = nHits();
  fwrite(&nh, sizeof(int), 1, fp);

  fwrite(&hits_[0], sizeof(Hit), nh, fp);
}

void Track::read_in(FILE *fp)
{
  fread(this, sizeof(Track), 1, fp);

  int nh = nHits();
  fread(&nh, sizeof(int), 1, fp);

  hits_.resize(nh);
  fread(&hits_[0], sizeof(Hit), nh, fp);
}
#endif
