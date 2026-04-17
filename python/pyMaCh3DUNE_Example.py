from pyMaCh3_DUNE import parameters
from pyMaCh3_DUNE import samples

from pathlib import Path
import os

if __name__=="__main__":
    # We can set up MaCh3 in a few ways!
    # The simplest is to just use the factory method
    
    
    config_file = Path("Configs")/"EventRates_Beam.yaml"
    if not config_file.is_file():
        raise FileExistsError(f"Cannot find {config_file}")
    
    # We then just set up all the parameter handlers in here!
    par_handler, sample_handlers = samples.MaCh3DuneFactory(str(config_file))
    