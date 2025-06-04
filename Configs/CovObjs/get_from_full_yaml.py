import yaml

def get_from_full_yaml(out_name,
                       det_id=None,
                       syst_type=None,
                       horn_current=None,
                       uniform_stepscale=False,
                       param_groups=None,
                       ):
    if det_id is None:
        det_id = [1]
    if syst_type is None:
        syst_type = ['Norm']
    if horn_current is None:
        horn_current = ['fhc','rhc']
    if param_groups is None:
        param_groups = ['Flux','Xsec','DetSys']
    with open("liban_grouped.yaml", 'r') as stream:
        try:
            f = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    new_yaml = {'Systematics': []}
    # b_systamatics = {1: [], 24: []}
    b_systamatics = []

    for s in f['Systematics']:
        # print("Systematic name: ", s['Systematic']['Names']['ParameterName'])
        accept = True
        if s['Systematic']['DetID'] not in det_id:
            accept = False
        # print(f'DetID: {s["Systematic"]["DetID"]} - {accept}')
        if s['Systematic']['Type'] not in syst_type:
            accept = False
        # print(f'Type: {s["Systematic"]["Type"]} - {accept}')
        if s['Systematic']['ParameterGroup'] not in param_groups:
            accept = False
        try:
            fhc_lower_bound = s['Systematic']['KinematicCuts'][1]['IsFHC'][0]
        except KeyError:
            fhc_lower_bound = 0
        if fhc_lower_bound >= 0:
            if 'fhc' not in horn_current:
                accept = False
        if fhc_lower_bound <= 0:
            if 'rhc' not in horn_current:
                accept = False
        # print(f'Horn current: {fhc_lower_bound} - {accept}')
        if accept:
            new_yaml['Systematics'].append(s)
            if 'b_' in s['Systematic']['Names']['ParameterName']:
                # detid = s['Systematic']['DetID'] 
                b_systamatics.append(s['Systematic']['Names']['ParameterName'])
            # else:
            #     new_yaml['Systematics'].append(s)
            #     if 'b' in s['Systematic']['Names']['ParameterName']:
            #         b_systamatics.append(s['Systematic']['Names']['ParameterName'])
    # Remove the FD systematics from the correlation
    # print(b_systamatics)
    for i,s in enumerate(new_yaml['Systematics']):
        if s['Systematic']['Names']['ParameterName'] in b_systamatics:
            corr = s['Systematic']['Correlations']
            corr = [c for c in corr if list(c.keys())[0] in b_systamatics]
            new_yaml['Systematics'][i]['Systematic']['Correlations'] = corr

        # Assume all norm systematics are flux systematics
        # if s['Systematic']['Type'] == 'Norm':
        #     new_yaml['Systematics'][i]['Systematic']['ParameterGroup'] = 'Flux'
        # elif s['Systematic']['Type'] == 'Spline':
        #     new_yaml['Systematics'][i]['Systematic']['ParameterGroup'] = 'Xsec'
        # elif s['Systematic']['Type'] == 'Functional':
        #     new_yaml['Systematics'][i]['Systematic']['ParameterGroup'] = 'DetSys'

        if uniform_stepscale:
            new_yaml['Systematics'][i]['Systematic']['StepScale']['MCMC'] = 1.0

    if len(new_yaml['Systematics']) > 0:
        with open(out_name, 'w') as stream:
            yaml.dump(new_yaml, stream)

if __name__ == '__main__':
    # get_from_full_yaml("liban_ND+FD_norm+spline_FHC+RHC.yaml",
    #                    det_id=[1,24,25],
    #                    syst_type=['Norm','Spline'],
    #                    horn_current=['fhc','rhc'])
    det_ids = [1,24,25]
    syst_types = ['Norm','Spline','Functional']
    param_groups = ['Flux','Xsec','DetSys']
    horn_currents = ['fhc','rhc']

    for det_id in det_ids:
        det_id_name = 'ND' if det_id == 1 else ('FD' if det_id == 24 else 'ND+FD')
        for param_group in param_groups:
            get_from_full_yaml(f"liban_covs/{param_group}_{det_id_name}.yaml",
                                det_id=[det_id],
                                syst_type=syst_types,
                                horn_current=horn_currents,
                                param_groups=[param_group])
            get_from_full_yaml(f"liban_covs/{param_group}_{det_id_name}_uniform.yaml",
                                det_id=[det_id],
                                syst_type=syst_types,
                                horn_current=horn_currents,
                                param_groups=[param_group],
                                uniform_stepscale=True)

    get_from_full_yaml('liban_covs/Flux_ND+FD.yaml',
                          det_id=[1,24],
                          syst_type=syst_types,
                          horn_current=['fhc','rhc'],
                          param_groups=['Flux'])
    get_from_full_yaml('liban_covs/Flux_ND+FD_uniform.yaml',
                       det_id=[1,24],
                          syst_type=syst_types,
                            horn_current=['fhc','rhc'],
                            param_groups=['Flux'],
                            uniform_stepscale=True)