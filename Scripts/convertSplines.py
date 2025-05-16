#!/usr/bin/python
#
# script to convert spline .root files into CUDA compatible classes

import os
import sys
import fnmatch
from ROOT import TFile, TIter, TKey, TSpline3

if (len(sys.argv) != 3):
    print "Useage: convertSplines.py /path/to/splines /path/to/output"
    quit(1)

### main 
data_loc = sys.argv[1]
output_loc = sys.argv[2] + "/"

# change these two lines to suit the files
match = "spl*numuoa*32*.root" # edit here
preprop = "double spl_"

os.system("mkdir -p spline_temp")

# make the base interface
base = open(output_loc + "splineInterface.h", "w")
base.write("#ifndef _splineInterface_h_\n")
base.write("#define _splineInterface_h_\n\n")
base.write("class splineInterface\n{\n")
base.write("public:\n")
base.write("virtual double fetchWeight(char *key, double frac)=0;\n};\n\n#endif")
base.close()

print "-----------------------------------------"
print " Converting TSpline3 objects to classes"
print "-----------------------------------------"

for f in os.listdir(data_loc):
    if (fnmatch.fnmatch(f, match)):
        print f
        os.system("mkdir -p spline_temp/" + f[:-5])
        
        file_o_functions = ""
        
        file1 = TFile(data_loc + "/" + f)
        #file1.ls()
        nextkey = TIter(file1.GetListOfKeys())
        
        # make init function
        content = "#include \"" + f[:-5] + ".h\"\n\n"
        content += "__host__ __device__ double " + str(f[:-5]) + "::fetchWeight(char *key, double frac)\n{\n"
        content += "if (key == "")\n{}\n"

        hList = ""

        for key in nextkey:
            if (str(key.GetClassName()) == "TSpline3"):
                
                #print str(key.GetClassName())
                spln = key.ReadObj()
                spl_name = str(spln.GetName())
                
                # generate splinename map entry for init function
                content += "else if (key == \"" + spl_name + "\") \n"
                content += "\treturn " + spl_name + "(frac);\n\tbreak;\n"

                spln.SaveAs(spl_name + ".cu")
                os.system("mv " + spl_name + ".cu spline_temp/" + f[:-5])
                
                # read the file back in
                txt = open("spline_temp/" + f[:-5] + "/" + spl_name + ".cu","r")
                funct = txt.read()
                funct = funct.replace(preprop, "__host__ __device__ " + preprop)
                
                funct = funct.replace("__host__ __device__ double ", "__host__ __device__ double " + f[:-5] + "::")

                hList += "__host__ __device__ double " + spl_name + "(double x);\n" 

                file_o_functions += funct + "\n"

                txt.close()
                # delete temp file
                #os.system("rm "+ spl_name + ".cpp")

        # finish switch function
        content += "default :\nreturn 1;\nbreak;\n"
        content += "}\n}\n\n"

        # save 
        outtxt = open(output_loc + f[:-5] + ".cu", "w")
        outtxt.write(content)
        outtxt.write(file_o_functions)
        outtxt.close()

        # make the header file
        header = open(output_loc + f[:-5] + ".h", "w")
        header.write("#ifndef _" + f[:-5] + "_h_\n")
        header.write("#define _" + f[:-5] + "_h_\n")
        header.write("\n#include \"splineInterface.h\"\n")
        header.write("\nclass " + f[:-5] + " : splineInterface\n{\npublic:\n")
        header.write("double fetchWeight(char *key, double frac);\n")
        header.write("protected:\n")
        header.write(hList)
        header.write("\n};\n\n#endif")
        header.close()
        
        os.system("rm -r spline_temp")
        
#os.system("mv " + spl_name + ".cpp spline_temp/" + f[:-5])
