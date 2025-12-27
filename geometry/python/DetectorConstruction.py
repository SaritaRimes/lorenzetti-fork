__all__ = [
           "DetectorConstruction", 
           ]

from typing import List

from GaugiKernel.constants import *
from GaugiKernel import Cpp, Logger
from GaugiKernel.macros import *

from prettytable import PrettyTable
from tqdm import tqdm

import os
import numpy as np
import collections

import ROOT

from geometry import PhysicalVolume, Plates
# import samplings
#from Tracking import createTracking
from ATLAS.ECAL import *
from ATLAS.TILE import *
from ATLAS.EMEC import *
from ATLAS.HEC  import *
from ATLAS.DeadMaterials import *
from ATLAS.Tracking      import *


basepath = os.environ['LORENZETTI_ATLAS_DATA_DIR']
vispath = f'{basepath}/vis.mac'

vis_begin = """
/vis/open OGL 600x600-0+0
/vis/viewer/set/autoRefresh false
/vis/verbose errors
/vis/drawVolume
/vis/viewer/set/viewpointVector 1 0 0
/vis/viewer/set/lightsVector 1 0 0
/vis/viewer/set/style wireframe
/vis/viewer/set/auxiliaryEdge true
/vis/viewer/set/lineSegmentsPerCircle 100
/vis/scene/add/trajectories smooth
/vis/modeling/trajectories/create/drawByCharge
/vis/modeling/trajectories/drawByCharge-0/default/setDrawStepPts true
/vis/modeling/trajectories/drawByCharge-0/default/setStepPtsSize 2
/vis/scene/endOfEventAction accumulate
/vis/geometry/set/visibility World 0 false
#/vis/geometry/set/visibility World 0 true
/vis/ogl/set/displayListLimit 10000000
"""
  
vis_command = """
/vis/geometry/set/colour {name} 0 {color}
/vis/geometry/set/colour {name}_Layer 0 {color}
/vis/geometry/set/colour {name}_Abso 0 {color}
/vis/geometry/set/colour {name}_Gap 0 {color}
/vis/geometry/set/visibility {name} 0 {visualization}
/vis/geometry/set/visibility {name}_Layer 0 {visualization}
/vis/geometry/set/visibility {name}_Abso 0 {visualization}
/vis/geometry/set/visibility {name}_Gap 0 {visualization}
"""

vis_end = """
/vis/viewer/set/autoRefresh true
/vis/verbose warnings
"""

def create_vis_mac(volumes, opath ):
  with open(opath, 'w') as f:
    f.write(vis_begin)
    for vol in volumes:
      f.write(vis_command.format(color=vol.Color, name=vol.name(), visualization='true' if vol.Visualization else 'false'))
    f.write(vis_end)

    

class DetectorConstruction( Cpp ):

  def __init__( self, 
                name              : str, 
                UseMagneticField  : bool=False, 
                CutOnPhi          : bool=False,
                vis_path          : str=vispath):

    Cpp.__init__(self, ROOT.DetectorConstruction(name) ) 
    
    self.setProperty( "UseMagneticField", UseMagneticField  )
    self.setProperty( "CutOnPhi"        , CutOnPhi          )


    samplings = []; volumes = []; tracking = []
    # Center
    #volumes.extend( getPixelBarrelCfg()   )
    samplings.extend( getLArBarrelCfg()   )
    samplings.extend( getTileBarrelCfg()  )
    volumes.extend( getDMVolumesCfg()      )
    # Right side (A)
    samplings.extend( getTileExtendedCfg()    )
    samplings.extend( getLArEMECCfg()         ) 
    samplings.extend( getHECCfg()             )
    volumes.extend( getCrackVolumesCfg()      )
    # Left side (B)
    samplings.extend( getTileExtendedCfg(left_side = True) )
    samplings.extend( getLArEMECCfg(left_side=True)        ) 
    samplings.extend( getHECCfg(left_side=True)            )    

    self.VisMac = vis_path
    self.__volumes = collections.OrderedDict({})

    # add all volumes from samplings
    for samp in volumes:
      pv = samp.volume(); self+=pv
  
  def __add__(self, pv):
    if pv.Name not in self.__volumes.keys():
      self.__volumes[pv.Name] = pv
    return self

  #
  # Build the core object
  #
  def compile(self):
    # Create all volumes inside of the detector
    for pv in tqdm( self.__volumes.values(), desc="Compiling...", ncols=70):
      self._core.AddVolume( pv.name(), pv.Plates, pv.AbsorberMaterial, pv.GapMaterial, 
                             # layer
                             pv.NofLayers, 
                             pv.AbsorberThickness, 
                             pv.GapThickness,
                             # dimensions
                             pv.RMin, pv.RMax, pv.ZSize, 
                             pv.X, pv.Y, pv.Z,
                             # production cuts
                             pv.Cuts.ElectronCut, 
                             pv.Cuts.PositronCut, 
                             pv.Cuts.GammaCut, 
                             pv.Cuts.PhotonCut
                             )
    #create_vis_mac(self.__volumes.values(), self.VisMac)


  def summary(self):

      samp_vol_names = []

      print('Display all calorimeter samplings...')

      t = PrettyTable(["Name", "Plates", "z",'Zmin','Zmax', "Rmin", 
                       "Rmax", "abso","gap", "deta", "dphi", "EtaMin", 
                       "EtaMax", "N_bins", "Container"])

      # Add all volumes that came from a sampling detector and has a sensitive parameter
      for samp in self.__volumes.values():
        pv = samp.volume(); sv = samp.sensitive(); samp_vol_names.append(pv.Name)
        t.add_row( [pv.Name,
                    Plates.tostring(pv.Plates),pv.ZSize,pv.ZMin,pv.ZMax,pv.RMin,pv.RMax,
                    pv.AbsorberMaterial,pv.GapMaterial,
                    round(sv.DeltaEta,4) ,
                    round(sv.DeltaPhi,4) ,
                    sv.EtaMin,sv.EtaMax, 
                    len(sv.EtaBins)*len(sv.PhiBins) if sv.DeltaEta else len(sv.ZBins)*len(sv.PhiBins)  ,
                    samp.CollectionKey
                  ])
      print(t)


      print('Display all tracking samplings...')

      t = PrettyTable(["Name", "Plates", "z",'Zmin','Zmax', "Rmin", 
                       "Rmax", "abso","gap", "dz", "dphi", "N_bins", "Container"])

      # Add all volumes that came from a sampling detector and has a sensitive parameter
      #for samp in self.samplings:
      #  pv = samp.volume(); sv = samp.sensitive(); samp_vol_names.append(pv.Name)
      #  t.add_row( [pv.Name,
      #              Plates.tostring(pv.Plates),pv.ZSize,pv.ZMin,pv.ZMax,pv.RMin,pv.RMax,
      #              pv.AbsorberMaterial,pv.GapMaterial,
      #              round(sv.DeltaEta,4) ,
      #              round(sv.DeltaPhi,4) ,
      #              sv.EtaMin,sv.EtaMax, 
      #              len(sv.EtaBins)*len(sv.PhiBins) if sv.DeltaEta else len(sv.ZBins)*len(sv.PhiBins)  ,
      #              samp.CollectionKey
      #            ])
      print(t)

      print('Display all non-sensitive volumes...')

      t = PrettyTable(["Name", "Plates", "z",'Zmin','Zmax', "Rmin", 
                       "Rmax", "abso","gap"])

      # Add ither volumes that not came from a sampling detector (extra volumes only)
      for key, pv in self.__volumes.items():
        if key not in samp_vol_names:
          t.add_row([pv.Name, Plates.tostring(pv.Plates),pv.ZSize, pv.ZMin, pv.ZMax, 
                     pv.RMin, pv.RMax, pv.AbsorberMaterial, pv.GapMaterial]) 
      print(t)














