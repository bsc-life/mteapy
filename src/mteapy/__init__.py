# Marks mteapy as a regular package (not a PEP 420 namespace package) --
# without this, Python merges in any other "mteapy" directory found earlier
# on sys.path (e.g. a sibling git checkout) alongside the real install.
