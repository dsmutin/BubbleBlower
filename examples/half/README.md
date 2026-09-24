# Half-strain assembly graph

Reads are the InSilicoSeq pair from `examples/half_strains` (samovar ISS).
metaSPAdes 4.2.0 from the existing conda env `vaegbin_env` writes
`work/metaspades`. The assembly graph is coloured by genome accession parsed
from the read ids, then BubbleBlower resolves bubbles.

`work/` is generated and is not part of the git tree. Run `run.py` after the
metaSPAdes directory contains `assembly_graph_with_scaffolds.gfa`.
