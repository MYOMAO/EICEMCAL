#pragma once

#include "GlobalVariables.C"
#include <g4main/PHG4Reco.h>

R__LOAD_LIBRARY(libg4testbench.so)

namespace G4WORLD
{
  double AddSpace = 100.;            // add this much space in cm around edge of detector
  string WorldMaterial = "G4_Galactic";   // default world material, use G4_Galactic for material scan
  string PhysicsList = "QGSP_BERT";  // for calorimeters use HP lists
}  // namespace G4WORLD

void WorldInit(PHG4Reco *g4Reco)
{
  g4Reco->SetWorldMaterial(G4WORLD::WorldMaterial);
  g4Reco->SetPhysicsList(G4WORLD::PhysicsList);
}


void WorldSize(PHG4Reco *g4Reco, double radius)
{
  double world_radius = radius;
  if (world_radius < BlackHoleGeometry::max_radius)
  {
    world_radius = BlackHoleGeometry::max_radius;
  }
  if (g4Reco->GetWorldSizeX() < (world_radius+100)) // add 1m
  {
    g4Reco->SetWorldSizeX(world_radius+100);
  }
// our world is a symmetric cylinder so the center is at 0/0/0, pick the largest of abs(min_z) || abs(max_z)
  double min_zval = std::min(BlackHoleGeometry::min_z,-((g4Reco->GetWorldSizeZ()-100)/2.));
  double max_zval = std::max(BlackHoleGeometry::max_z,(g4Reco->GetWorldSizeZ()-100)/2.);
  double final_zval = std::max(fabs(min_zval),fabs(max_zval))+100; // add 1m
  if (g4Reco->GetWorldSizeZ() < 2*(final_zval))
  {
    g4Reco->SetWorldSizeZ( 2*(final_zval));
  }
  return;
}
