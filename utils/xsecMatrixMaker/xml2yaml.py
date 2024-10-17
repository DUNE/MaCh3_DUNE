#!/usr/bin/env python3

import sys
import yaml
import xml.etree.ElementTree as ET

def convert_parameter(param:"dict[str, str]") -> dict:
    yaml_param = {}

    yaml_param['Names'] = {}
    yaml_param['Names']['ParameterName'] = param['name']
    yaml_param['Names']['FancyName'] = param['name']

    yaml_param['ParameterValues'] = {}
    yaml_param['ParameterValues']['PreFitValue'] = float(param['prior']) #TODO: check which is which between prior and nom
    yaml_param['ParameterValues']['Generated'] = float(param['nom'])

    yaml_param['ParameterBounds'] = [float(param['lb']), float(param['ub'])]

    yaml_param['StepScale'] = {}
    yaml_param['StepScale']['MCMC'] = float(param['stepscale'])

    yaml_param['DetID'] = int(param['detid'])
    yaml_param['Error'] = float(param['error'])
    yaml_param['Type'] = param['type'].capitalize()
    yaml_param['FlatPrior'] = False

    if 'Spline' in yaml_param['Type']:
        yaml_param['SplineInformation'] = {}
        if 'fd_spline_name' in param:
            yaml_param['SplineInformation']['SplineName'] = param['fd_spline_name']
        if 'fd_mode' in param:
            yaml_param['SplineInformation']['FDMode'] = list(map(int, param['fd_mode'][0]['value'].split()))

    if 'kinematic_cut' in param:
        yaml_param['KinematicCuts'] = []
        if not isinstance(param['kinematic_cut'], list):
            param['kinematic_cut'] = [param['kinematic_cut']]
        for cut in param['kinematic_cut']:
            yaml_param['KinematicCuts'].append({cut['var']: list(map(float, cut['value'].split()))})

    if 'correlation' in param:
        yaml_param['Correlations'] = []
        if not isinstance(param['correlation'], list):
            param['correlation'] = [param['correlation']]
        for corr in param['correlation']:
            yaml_param['Correlations'].append({corr['par']: float(corr['value'])})

    return {'Systematic': yaml_param}



def xml2yaml(ifile:str, ofile:str) -> None:
    tree = ET.parse(ifile)
    root = tree.getroot()

    assert(root.tag == 'data')

    systematics = []

    for child in root:
        attrib = child.attrib
        for subelt in child:
            sub_attrib = subelt.attrib 
            sub_attrib.update({'value': subelt.text})

            if subelt.tag not in attrib:
                attrib[subelt.tag] = []
                
            attrib[subelt.tag].append(sub_attrib)

        syst = convert_parameter(attrib)
        systematics.append(syst)
    
    with open(ofile, 'w') as f:
        yaml.safe_dump({'Systematics': systematics}, f)

    

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} ifile.xml ofile.yaml")
        sys.exit(1)
    xml2yaml(sys.argv[1], sys.argv[2])
