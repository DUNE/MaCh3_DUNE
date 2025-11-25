import uproot
import numpy as np
import sys

fn = sys.argv[1]
with uproot.open(fn) as f:
    print("Keys:", f.keys())
    for key in ("correlation_matrix","Correlation","covariance_matrix"):
        if key in f:
            print("Found", key)
            arr = f[key].array() if hasattr(f[key], "array") else None
            # uproot may expose TMatrix as raw; try extracting by reading as numpy via tolist
            try:
                mat = f[key].to_numpy()  # uproot 5-ish; might fail on some versions
            except Exception:
                try:
                    mat = np.array(f[key].read())
                except Exception:
                    # fallback: read element-wise if it's saved as a TH2 or array
                    print("Could not auto-read matrix for", key)
                    mat = None
            print("matrix type:", type(mat))
            if isinstance(mat, np.ndarray):
                print("shape:", mat.shape)
                print("min/max:", np.nanmin(mat), np.nanmax(mat))
                print("any NaN? ", np.isnan(mat).any())
                print("any Inf? ", np.isinf(mat).any())
            else:
                print("Manual check required for", key)
