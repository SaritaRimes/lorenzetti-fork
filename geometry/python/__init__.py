
__all__ = ['flatten']



def flatten(l):
    o = []
    for item in l:
        if type(item) is list:
            o.extend(item)
        else:
            o.append(item)
    return o


from . import PhysicalVolume
__all__.extend(PhysicalVolume.__all__)
from .PhysicalVolume import *

from . import SensitiveDetector
__all__.extend(SensitiveDetector.__all__)
from .SensitiveDetector import *

from . import DetectorConstruction
__all__.extend(DetectorConstruction.__all__)
from .DetectorConstruction import *

from .detectors import *
__all__.extend(detectors.__all__)
from .detectors import *

