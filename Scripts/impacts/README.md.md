# Impact Plots

See this talk for full description of impact plots: https://indico.fnal.gov/event/71671/timetable/?view=standard#81-validity-of-systematics-mod

The code originally used to make the plots in that code was a little hacky. The code found in this repository is less hacky and more generic, but makes plots that are not as pretty. Furthermore, the feature of plotting pulls in different windows of the parameter phase space has not been re-implemented here yet. These issues can be solved in the future.

### Make Impacts

To create an impacts plots, first use the `make_impacts.py` script to create a `impacts.json`. Usage can be found with `make_impacts.py --help`:
```
(MaCh3) [scarf1534@cn1052 impacts]$ python3 make_impacts.py --help
usage: make_impacts.py [-h] [--metric-param METRIC_PARAM] [--save-path SAVE_PATH] [--burn-in BURN_IN] chain_path metric_type

positional arguments:
  chain_path
  metric_type

options:
  -h, --help            show this help message and exit
  --metric-param, -p METRIC_PARAM
  --save-path SAVE_PATH
  --burn-in BURN_IN
```

Minimally, you need to specify the MCMC chain ROOT file and the name of the metric to evaluate impact on. For most metrics, you also have to specify a parameter for the metric. For example:
```
python3 make_impacts.py my_chain.root Mean --metric-param delta_cp
```
which reads `my_chain.root`, and evaluates the impact of the Median value of the 1D $\delta_{CP}$ posterior.

All avaliable metrics are:
- Mean
- HPD
- Median
- StdDev
- IntervalWidth68
- IntervalMid68
- IntervalWidth95
- IntervalMid95
- MassOrderingBF
- MassOrderingBFSigma

All except MassOrderingBF and MassOrderingBFSigma require the metric parameter to be specified with `metric-param`.

> [!CAUTION]
> A reweighting factor is currently hard-coded into MassOrderingBF and MassOrderingBFSigma which assumes a penalty of exp(10) was applied to the likelihood for NO steps.

> [!CAUTION]
> The metrics that rely on binning, like HPD and IntervalWidth/Mid can create weird impacts which stem from binning dependence. For the time being, I'd recommend unbinned metrics like Mean and StdDev until more studies have been done.

You may want to set a burn in for the chain as well which you can do with the `--burn-in` option. For example:
```
python3 make_impacts.py my_chain.root Mean --metric-param delta_cp --burn-in 100000
```

### Make Pulls

We now use the `make_pulls.py` script to create `pulls.json`. The script just requires the MCMC chain ROOT file:
```
python3 make_pulls.py my_chain.root
```

### Plotting

Now we bring `impacts.json` and `pulls.json` together to produce `impacts.pdf`
```
python plot_impacts.py impacts.json pulls.json
```

Example plot:
![impact plot](image.png)