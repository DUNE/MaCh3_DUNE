'''
Very simple script that just loops over files
'''

import uproot
import glob
from rich import print

def check_corrupt(file: str):
    with uproot.open(file) as f:
        print(f"[blue]{file}[/blue] [bold green]FINE")

def main(file_folder: str):

    fail = 0
    
    for file in glob.glob(f"{file_folder}/*.root"):
        try:
            check_corrupt(file)
        except Exception:
            fail += 1 
            print(f"[yellow]{file}[/yellow] [bold red]CORRUPT!")
        
    if fail>0:
        raise Exception(f"[yellow]Failed check on [/yellow] [bold red]{fail} files")
        
if __name__=="__main__":
    import sys

    for dir in sys.argv[1:]:
        main(dir)