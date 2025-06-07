import yaml
import argparse

def get_from_full_yaml(input_file, out_name,
                       det_id=None,
                       syst_type=None,
                       horn_current=None,
                       uniform_stepscale=False,
                       param_groups=None,
                       ):
    if det_id is None:
        det_id = [["ND"]]
    if syst_type is None:
        syst_type = ['Norm']
    if horn_current is None:
        horn_current = ['fhc','rhc']
    if param_groups is None:
        param_groups = ['Flux','Xsec','DetSysts']
    with open(input_file, 'r') as stream:
        try:
            f = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    new_yaml = {'Systematics': []}
    b_systamatics = []

    print("Looping over systematics")
    for s in f['Systematics']:
        if set(s['Systematic']['DetID']) not in det_id:
            continue
        if s['Systematic']['Type'] not in syst_type:
            continue
        if s['Systematic']['ParameterGroup'] not in param_groups:
            continue
        try:
            fhc_lower_bound = s['Systematic']['KinematicCuts'][1]['IsFHC'][0]
        except KeyError:
            fhc_lower_bound = 0
        if fhc_lower_bound >= 0:
            if 'fhc' not in horn_current:
                continue
        if fhc_lower_bound <= 0:
            if 'rhc' not in horn_current:
                continue
        new_yaml['Systematics'].append(s)
        if 'b_' in s['Systematic']['Names']['ParameterName']:
            b_systamatics.append(s['Systematic']['Names']['ParameterName'])
    print("Filling new yaml")
    for i,s in enumerate(new_yaml['Systematics']):
        if s['Systematic']['Names']['ParameterName'] in b_systamatics:
            corr = s['Systematic']['Correlations']
            corr = [c for c in corr if list(c.keys())[0] in b_systamatics]
            new_yaml['Systematics'][i]['Systematic']['Correlations'] = corr

        if uniform_stepscale:
            new_yaml['Systematics'][i]['Systematic']['StepScale']['MCMC'] = 1.0

    if len(new_yaml['Systematics']) > 0:
        with open(out_name, 'w') as stream:
            yaml.dump(new_yaml, stream)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Split full yaml into smaller yaml files')
    parser.add_argument('input', type=str, help='Input yaml file')
    parser.add_argument('output_folder', type=str, help='Output folder')
    
    args = parser.parse_args()
    input_file = args.input
    output_folder = args.output_folder
    
    det_ids = [set(['ND']), set(['FD']), set(['ND','FD'])]
    det_ids_names = ['ND', 'FD', 'ND+FD']
    syst_types = ['Norm','Spline','Functional']
    param_groups = ['Flux','Xsec','DetSysts']
    horn_currents = ['fhc','rhc']

    for i,det_id in enumerate(det_ids):
        det_id_name = det_ids_names[i]
        for param_group in param_groups:
            get_from_full_yaml(input_file, f"{output_folder}/{param_group}_{det_id_name}.yaml",
                                det_id=[det_id],
                                syst_type=syst_types,
                                horn_current=horn_currents,
                                param_groups=[param_group])
            get_from_full_yaml(input_file, f"{output_folder}/{param_group}_{det_id_name}_uniform.yaml",
                                det_id=[det_id],
                                syst_type=syst_types,
                                horn_current=horn_currents,
                                param_groups=[param_group],
                                uniform_stepscale=True)

    # Combine ND and FD for flux
    get_from_full_yaml(input_file, f'{output_folder}/Flux_ND+FD.yaml',
                          det_id=[set(['ND']),set(['FD'])],
                          syst_type=syst_types,
                          horn_current=['fhc','rhc'],
                          param_groups=['Flux'])
    get_from_full_yaml(input_file, f'{output_folder}/Flux_ND+FD_uniform.yaml',
                       det_id=[set(['ND']),set(['FD'])],
                          syst_type=syst_types,
                            horn_current=['fhc','rhc'],
                            param_groups=['Flux'],
                            uniform_stepscale=True)
