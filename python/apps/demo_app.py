from _pyMaCh3 import MaCh3Instance
import argparse

# HW : Really really simple python app

if __name__=="__main__":
    parser = argparse.ArgumentParser(usage="python demo_app -c <config_name>.yaml")
    parser.add_argument("-c", "--config", help="TOML config file")
    
    args = parser.parse_args()
    
    mach3 = MaCh3Instance(args.config)
    
    print(f"Parameter values are : {mach3.get_parameter_values}")
    
    
    parameter_values[0] += 1
    print(f"Likelihood after shift : {mach3.propose_step(parameter_values)}")
    

    