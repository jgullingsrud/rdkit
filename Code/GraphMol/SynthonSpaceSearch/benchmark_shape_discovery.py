"""Compare native shape-SMILES discovery with an embedding-free callback.

Use a shape-free binary .spc database, and run each mode in a separate process.
The checksum covers sorted (SMILES, requested conformer count) pairs.
"""

import argparse
import hashlib
import json
import time

from rdkit import Chem, rdBase
from rdkit.Chem import rdEnumerateStereoisomers, rdSynthonSpaceSearch


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("database")
  parser.add_argument("--legacy", action="store_true")
  parser.add_argument("--threads", type=int, default=8)
  parser.add_argument("--max-synthon-atoms", type=int, default=50)
  args = parser.parse_args()
  if args.threads < 1 or args.max_synthon_atoms < 0:
    parser.error("threads must be positive and the atom limit nonnegative")
  params = rdSynthonSpaceSearch.ShapeBuildParams()
  params.maxSynthonAtoms = args.max_synthon_atoms
  params.interimWrites = 0
  options = rdEnumerateStereoisomers.StereoEnumerationOptions()
  options.randomSeed = params.stereoEnumOpts.randomSeed
  options.tryEmbedding = False
  params.stereoEnumOpts = options
  started = time.perf_counter()
  space = rdSynthonSpaceSearch.SynthonSpace()
  with rdBase.BlockLogs():
    space.ReadDBFile(args.database, args.threads)
    loaded = time.perf_counter()
    if space.GetNumSynthonsWithShapes():
      parser.error("input must be shape-free")
    if args.legacy:
      requests = set()

      def record(smiles, num_confs):
        requests.add(smiles)
        return Chem.MolFromSmiles(smiles)

      params.numThreads = 1
      params.setUserConformerGenerator(record)
      space.BuildSynthonShapes(params)
      smiles = sorted(requests)
    else:
      params.numThreads = args.threads
      smiles = space.GetSynthonShapeSmiles(params)
  finished = time.perf_counter()
  digest = hashlib.sha256(b"[")
  for index, smile in enumerate(smiles):
    if index:
      digest.update(b", ")
    digest.update(json.dumps([smile, params.numConfs]).encode())
  digest.update(b"]")
  print(json.dumps({
    "mode": "legacy" if args.legacy else "native",
    "threads": params.numThreads,
    "synthons": space.GetNumSynthons(),
    "requests": len(smiles),
    "load_s": loaded - started,
    "discovery_s": finished - loaded,
    "total_s": finished - started,
    "sha256": digest.hexdigest(),
  }, indent=2))


if __name__ == "__main__":
  main()