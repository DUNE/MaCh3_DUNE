
from pathlib import Path
import os
from typing import TypedDict, Literal

from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from pyMaCh3_DUNE import parameters
from pyMaCh3_DUNE import samples


########### Sample setup methods ###############
SampleType= Literal['beam_fd', 'beam_nd', 'atm', 'beam_nd_gar']

class SampleInput(TypedDict):
    path: Path
    type: SampleType

def make_sample(sample_path: Path, sample_type: SampleType, parameter_handler: parameters.ParameterHandlerGeneric):
    '''
    Makes a sample handler object
    '''
    
    sample_dict: dict[SampleType, samples.SampleHandlerBase] = {"beam_fd": samples.SampleHandlerBeamFD,
                                                                "beam_nd": samples.SampleHandlerBeamND,
                                                                "atm": samples.SampleHandlerAtm,
                                                                "beam_nd_gar":  samples.SampleHandlerBeamNDGAr}

    if not sample_path.is_file():
        raise FileNotFoundError(f"Cannot find {sample_path}")
    
    sample_handler = sample_dict.get(sample_type)
    
    if (sample_handler) is None:
        raise KeyError(f"Cannot find '{sample_type}' in samples! ('{list(sample_dict.keys())}')")
    
    return sample_handler(str(sample_path), parameter_handler)

def _check_files_exist(files: list[Path]):
    not_found = []
    for p in files:
        if p.is_file:
            continue
        not_found.append(p)
    
    if not_found:
        raise FileNotFoundError(f"Cannot find config files: {not_found}")

def setup_parameters(mach3_path: Path)->parameters.ParameterHandlerGeneric:
    '''
    Sets up the parameter handler
    '''
    config_dir= mach3_path/"Configs"/"CovObjs"
    tdr_dir = config_dir/"tdr_covs"
    
    # Configs
    parameter_files = [
        tdr_dir/"Flux_ND+FD.yaml",
        tdr_dir/"Xsec_ND+FD.yaml",
        tdr_dir/"DetSysts_FD.yaml",
        config_dir/"OscCov_PDG2021_v2.yaml"
    ]
           
    _check_files_exist(parameter_files)
    
    return parameters.ParameterHandlerGeneric([str(p) for p in  parameter_files])

def setup_samples(mach3_path: Path, parameters: parameters.ParameterHandlerGeneric)->list[samples.SampleHandlerBase]:
    config_path = mach3_path/"Configs"/"Samples"
    sample_configs: list[SampleInput] = [
        {"path": config_path/"SampleHandler_FD_Beam.yaml",
         "type": "beam_fd"}]
    
    return [make_sample(s.get('path'), s.get('type'), parameters) for s in sample_configs]    
    
def plot_mc_hist(sample_handlers: list[samples.SampleHandlerBase], output_file: Path):
    '''
    Plots the MC histograms    
    '''

    output_file.parent.mkdir(exist_ok=True, parents=True)
    
    with PdfPages(str(output_file)) as pdf:
        for handler in sample_handlers:
            handler.reweight()
            for n in range(handler.get_n_samples()):
                mc_pred, bin_edges_x, bin_edges_y = handler.get_mc_hist(n)
                            
                # 1D hist ()

                if not len(bin_edges_y):
                    plt.hist(bin_edges_x[:-1], bins=bin_edges_x, weights=mc_pred, histtype="step")
                else:
                    plt.hist2d(bin_edges_x[:-1], bin_edges_y[:-1], [bin_edges_x, bin_edges_y], weights=mc_pred)
                pdf.savefig()
                plt.close()
    
if __name__=="__main__":
    # We can set up MaCh3 in a few ways!
    # The simplest is to just use the factory method
    mach3_dir = os.environ.get("MACH3")
    if mach3_dir is None:
        raise OSError("Error 'MACH3' not set. Go to your top level MaCh3 directory and `export MACH3=$(pwd)`")

    mach3_path = Path(mach3_dir)

    parameter_handler = setup_parameters(mach3_path)
    sample_handlers = setup_samples(mach3_path, parameter_handler)
    plot_mc_hist(sample_handlers, mach3_path/"mc_hist_example.pdf")