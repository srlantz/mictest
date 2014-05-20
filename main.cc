/*
  g++ -std=c++11 -o main main.cc Track.cc Hit.cc Matrix.cc KalmanUtils.cc Propagation.cc Simulation.cc buildtest.cc fittest.cc -I. `root-config --libs --cflags`
  icc -std=gnu++0x -O3 -openmp -o main main.cc Track.cc Hit.cc Matrix.cc KalmanUtils.cc Propagation.cc Simulation.cc buildtest.cc fittest.cc -I. `root-config --libs --cflags`
*/

#include <iostream>

#include "MatriplexCommon.h"

#include "fittest.h"
#include "buildtest.h"

int main()
{
  bool saveTree = true;

  double t0, tmp, tsm;

  t0 = dtime();
  runFittingTest(saveTree,5000);
  tsm = dtime() - t0;

  t0 = dtime();
  runFittingTestPlex(saveTree, 5000);
  tmp = dtime() - t0;

  printf("SMatrix = %.3f   Matriplex = %.3f   ---   SM/MP = %.3f\n", tsm, tmp, tsm / tmp);

  // runBuildingTest(saveTree,10);

  return 0;
}
