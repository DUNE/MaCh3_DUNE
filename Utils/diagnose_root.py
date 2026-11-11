#!/usr/bin/env python3
"""
diagnose_root.py -- dump the structure of a posterior-predictive ROOT file
so we can see exactly how it's laid out, rather than guessing.

Usage: python diagnose_root.py file.root
"""
import sys
import uproot

def walk(node, prefix="", depth=0, max_depth=5, max_keys_per_level=25):
    if depth > max_depth:
        return
    try:
        keys = node.keys(recursive=False)
    except Exception as exc:
        print("  " * depth + f"[error listing keys: {exc}]")
        return

    shown = keys[:max_keys_per_level]
    for k in shown:
        ks = k.split(";")[0]
        try:
            obj = node[k]
            is_hist = hasattr(obj, "axes")
            tag = "[HIST]" if is_hist else "[DIR]"
        except Exception as exc:
            tag = f"[error: {exc}]"
            obj = None
        print("  " * depth + f"{tag} {ks}")
        if obj is not None and not hasattr(obj, "axes"):
            walk(obj, depth=depth + 1, max_depth=max_depth,
                 max_keys_per_level=max_keys_per_level)
    if len(keys) > max_keys_per_level:
        print("  " * depth + f"... ({len(keys) - max_keys_per_level} more keys not shown)")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python diagnose_root.py file.root")
        sys.exit(1)
    f = uproot.open(sys.argv[1])
    print(f"=== {sys.argv[1]} ===")
    walk(f, max_depth=5, max_keys_per_level=25)