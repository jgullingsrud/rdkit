# Embedding-Free Shape Discovery

`SynthonSpace.GetSynthonShapeSmiles(params)` returns sorted, unique SMILES for
shape preparation without generating conformers or modifying the space. It
includes the stereo variants and unenumerated fallback molecules that an
embedding-free `BuildSynthonShapes` conformer callback would see.

```python
from rdkit.Chem import rdSynthonSpaceSearch

space = rdSynthonSpaceSearch.SynthonSpace()
space.ReadDBFile("REAL.spc", 8)
params = rdSynthonSpaceSearch.ShapeBuildParams()
params.numThreads = 8
params.maxSynthonAtoms = 50
smiles = space.GetSynthonShapeSmiles(params)
print(len(smiles))
```

The implementation uses native workers and thread-local result sets. It avoids
Python callbacks, SMILES reparsing, conformer matching, and shape-construction
retries. Exact assembly contexts (target role plus all partner synthon pointers)
are processed once when the stereo seed is fixed. This does not merge merely
similar structures or discard specified stereo. Both discovery and shape
building use the same stereoisomer limit and sample construction/trimming code.

`maxSynthonAtoms` and existing shapes are respected. `stereoEnumOpts` controls
enumeration, except that `tryEmbedding` is always disabled. Use a fixed
`stereoEnumOpts.randomSeed` for results independent of worker count. The conformer
generator is never called; conformer counts and embedding settings are not used.
Worker exceptions propagate to the caller. Cancellation raises `KeyboardInterrupt`
in Python without returning a partial result. Concurrent mutation of the space
is not supported. This is a batch API, not a streaming replacement for callbacks.

## Benchmark

Use a shape-free native binary database and run each mode in a fresh process:

```sh
python Code/GraphMol/SynthonSpaceSearch/benchmark_shape_discovery.py REAL.spc --legacy
python Code/GraphMol/SynthonSpaceSearch/benchmark_shape_discovery.py REAL.spc --threads 8
```

Compare both request counts and checksums. The benchmark hashes sorted
`(SMILES, numConfs)` pairs, not just the number of requests.

Measured on a 14-core ARM macOS machine, using RDKit 2026.03.6 and a release-mode
build of this change:

| Input and Operation | Elapsed | Unique Requests |
| --- | ---: | ---: |
| REAL 2025-09 compressed input, original loading + callback discovery | 2,888.52 s | 3,203,796 |
| Same input, complete updated counting command, 8 threads | 232.67 s | 3,203,796 |
| REAL 2025-09 native binary, new API benchmark, 8 threads | 226.58 s | 3,203,796 |

The complete compressed-input counting command was **12.4x faster**. Its peak RSS
was about 15.0 GB. The baseline harness peaked at 12.2 GB, including checksum
construction; these peak-memory measurements include different output handling.
These are single full-size runs, not repeated isolated benchmark averages. Some
experiments overlapped the baseline run; repeat measurements on deployment
hardware for capacity planning. The original timing excludes interpreter startup
and final checksum construction; the updated command timing includes process
startup and teardown.

The input contains 552,083 distinct synthons in 1,015 reactions. Exact-context
deduplication reduced 7,281,671 candidates to 6,400,340. Complete baseline,
prototype, and integrated-API request checksums matched:

```text
ac9adaa017ed7f4f20e557d2db5f9f2781496c7365b7cafd571a1282e461bd85
```

On the 50,117-synthon fixture, native API discovery took 8.14 / 1.91 / 1.57 seconds
with 1 / 8 / 14 threads. All produced the same 83,182 requests. These timings
exclude input loading.

Discarded experiments included skipping fallback based only on stereoisomer
count, which lost 295 fixture requests due to specified ring stereo, and in-place
sample assembly, which did not show a useful improvement. Shape-building retry
semantics are unchanged.

## Reusing Preparation During Assembly

`PrepareSynthonShapes(params, callback)` emits `(smiles, record)` for each
eligible synthon. The opaque byte record contains the already-trimmed molecules,
their atom/bond properties, precomputed request SMILES, and stereo/fallback
groups in the existing context-preference order. Records are prepared in bounded
batches using native threads; the callback runs on the calling thread.

`BuildSynthonShapesFromPreparation(params, supplier, conformerSupplier=None)`
consumes those records after the original binary database has been reloaded.
`supplier()` returns the next record or `None` at EOF. Conformers can come from
the usual user generator, or from the optional `(smiles, count)` supplier of
uncompressed RDKit molecule-pickle bytes (`None` means failure/missing).
Binary molecules are decoded outside the Python GIL. Prepared templates are
decoded only when conformers are available. Shape extraction, matching, pruning,
and the existing context/stereo fallback policy are shared with ordinary builds.

The records are versioned and checked against the RDKit version and preparation
settings. The caller must bind them to the exact source database using a checksum;
the floe does this and checks transferred record counts and byte totals. Record
delivery order does not select a context; preference order is stored inside each
synthon record. Duplicate synthon records are rejected. Exceptions from record
suppliers and worker-side Python conformer suppliers propagate to the caller.
The prepared API does not currently emit intermediate database snapshots or a
progress bar. The ordinary shape-building API retains those features.

The floe defaults to this path when supported and discovery is unbounded, spooling
compressed records on disk through its manifest stream. Older RDKit builds,
bounded discovery, and `reuse_discovery=False` retain the previous path.

Paired eight-thread assembly-and-write measurements, excluding loading and
embedding, with two alternating-order repetitions:

| Binary Input / Cache | Legacy | Prepared | Preparation Artifact |
| --- | ---: | ---: | ---: |
| 50,117 synthons / 24 successful requests | 9.17 s | 2.80 s | 46.3 MB |
| 2,000 synthons / all 3,400 requests computed | 0.53 s | 0.48 s | 1.86 MB |

Both used two requested conformers. The fully computed cache had 3,352 successes
and 48 failures, producing 1,999 shaped synthons. All final database checksums
matched their corresponding legacy builds, including tests with alternate
reaction failures. Preparation took 4.16 s and 0.177 s respectively. These
results do not establish a full-REAL assembly speedup; storage/transfer overhead
and high conformer coverage can substantially reduce the benefit.

Cold-start testing also exposed an existing race in GaussianShape's feature
pattern initialization: its pointer was published before population finished.
The local build includes a thread-safe static initializer and a concurrent-first-
use regression test. This changes no feature definitions or selection rules.

The tested assembly overlay is `/tmp/rdkit-synthon-perf/assembly-lazy`. Set
`RDKIT_SYNTHON_OVERLAY` and `DYLD_LIBRARY_PATH` to that directory rather than
`native-final` below when testing preparation reuse. It also contains the fixed
GaussianShape library. The installed wheel remains unchanged.

## Distributed Shape Construction

Three additional operations separate shape computation from database mutation:

* `GetSynthonShapePreparationInfo(params, record)` returns the target synthon key
	and its unique conformer dependencies without reconstructing molecules.
* `BuildSynthonShapeFromPreparation(params, record, conformerSupplier)` is a
	module-level function requiring no `SynthonSpace`. It returns the target key
	and serialized `SynthonShapeInput` bytes. Empty bytes mean no context succeeded.
	The supplier returns molecule-pickle bytes or `None` for unavailable conformers.
* `space.SetSynthonShapes(entries, numConformers)` attaches `(key, shapeBytes)`
	entries without calculating shapes. Unknown/repeated/already-shaped targets
	are rejected within a batch, and all payloads are parsed before applying it.

Independent workers use the same shape-building implementation as local assembly.
Each record contains the whole context preference sequence for one synthon, so
completion order does not change selection. The final writer only loads the
source database, attaches results, writes, and validates the ordinary `.spc` file.

The floe routes validated conformer results and preparation into deterministic
self-contained ZIP bundles. Only required conformers are included, deduplicated
within each job. Per-job manifests bind source identity, configuration, targets,
and preparation checksums. Identical worker retries are idempotent; missing jobs,
conflicting retries, unexpected targets, and provenance mismatches prevent output.
Failed/missing conformers remain nonfatal and retain the ordinary fallback rules.

The router waits for conformer generation to finish. It keeps payloads on disk,
then groups work using `shape_batch_size` and a soft `shape_batch_mb` target.
One synthon is indivisible. Conformers needed by different jobs are retransmitted;
routing, network traffic, and final serialization remain scaling limits. Shape
workers have no synthon-library input. `max_shape_workers` controls distributed
concurrency, and `shape_cpus` controls native work within each worker.

The 2,000-synthon binary sample was built in 63 jobs on two spawned processes
with two threads each. The final database was byte-identical to legacy assembly,
including 48 failed conformer requests and 1,999 successfully shaped synthons.
The local multi-process run took 1.23 s including routing, startup, attachment,
and writing, versus 0.75 s for ordinary four-thread assembly. This is a correctness
check, not a cloud speedup claim; full-REAL distributed throughput is unmeasured.

The tested overlay for these APIs is `/tmp/rdkit-synthon-perf/distributed-shapes`.
Use it for both `RDKIT_SYNTHON_OVERLAY` and the leading `DYLD_LIBRARY_PATH` entry.
The existing `overlay/sitecustomize.py` remains the Python import hook. No installed
wheel was modified. Deployment requires a rebuilt RDKit package providing these
APIs; unsupported builds, bounded runs, and `parallel_shapes=False` use local assembly.

## Local Development Build

The benchmark and runtime-test results above were obtained on `Release_2026_03_6`
with rebuilt synthon and GaussianShape libraries and the Boost.Python wrapper.
This was not a full RDKit build, and the installed wheel was not modified.

The implementation branch has since been rebased onto upstream `master` at
`cbfb37abddcd5b5feeac97d53530ae6be83cac0d`. The rebase preserves master's mutable
`trimSampleMol` signature, ring handling, and other upstream fixes. The rebased
native implementation, Boost.Python wrapper, and native tests were compile-checked;
the Python test sources were also compiled. Runtime tests and benchmarks have not
yet been rerun against a complete master-based build, and nanobind support for
these APIs has not been added. A new build must use RDKit's standard CMake
configuration, including threading and Boost serialization.

The temporary overlays and `master-forwardport.patch` mentioned above are
historical artifacts, not builds or patches of the rebased branch. In particular,
the release-based libraries must not be used to validate master's changed ABI.

To reproduce the earlier release-based experiment in this workspace, run from
`tl-synthon-searching`:

```sh
export PYTHONPATH=/tmp/rdkit-synthon-perf/overlay
export RDKIT_SYNTHON_OVERLAY=/tmp/rdkit-synthon-perf/native-final
export DYLD_LIBRARY_PATH="$RDKIT_SYNTHON_OVERLAY:$PWD/.venv/lib/python3.11/site-packages/rdkit/.dylibs"
.venv/bin/python scripts/count_discovery_smiles.py 2025-09_REAL_synthons.spc.zst
```

These temporary build artifacts are machine-specific. Without a build providing
the new API, the counting script falls back to its original callback path.

Native shape-search tests, Python wrapper coverage, and application shape tests
exercise the new API. A full RDKit test suite remains a CI requirement before
upstream submission. The change and this experiment used GitHub Copilot; disclose
that assistance in the PR description per RDKit's contribution guidelines.