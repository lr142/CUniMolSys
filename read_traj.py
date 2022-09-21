import sys
sys.path.append("../")
sys.path.append("D:/_Nutbox/UniversalMolecularSystem")
sys.path.append("/mnt/d/_Nutbox/UniversalMolecularSystem")
sys.path.append("/Users/me/Projects/UniversalMolecularSystem")

from UniversalMolecularSystem import *
from Utility import *
from LAMMPSDATAFile import *
from LAMMPSDUMPFile import *


def Main():
    PATH = "/Users/me/Projects/CUniMolSys/0402_PolyDADMAC_400K"
    ms = MolecularSystem()
    ms.Read(LAMMPSDATAFile(),os.path.join(PATH,"system.data"))
    ms.Summary()
    ms.ReadTrajectory(LAMMPSDUMPFile(),os.path.join(PATH,"system.lammpstrj"),maxFrames=99999,flushSameTimestep=True,
                      timestep_in_fs=1.0,every=1,max_workers=1)
    ms.ReadTrajectory(LAMMPSDUMPFile(), os.path.join(PATH, "system.lammpstrj.2"), maxFrames=999999, flushSameTimestep=True,
                      timestep_in_fs=1.0,every=1,max_workers=1)


Main()
