# denovo_OLC

A standalone release of the [LFR pipeline](https://github.com/Complete-Genomics/LFR_Pipeline) denovo module: an in-process Overlap-Layout-Consensus (OLC) assembler purpose-built for **single-end 600bp (SE600) linked-read data**, reconstructing one fragment per UMI/barcode at millions-of-molecules scale. Designed as a drop-in alternative to per-UMI de Bruijn graph assembly (e.g. forking a MEGAHIT subprocess per barcode), and engineered specifically around the noise and read-structure characteristics that SE600 data — as opposed to short paired-end reads — actually presents.

## Why SE600 needs a purpose-built assembler

Complete LFR (cLFR) and similar barcode-partitioned sequencing technologies produce, per physical DNA molecule, a small cluster of short reads sharing one barcode. Most existing per-barcode assembly tooling is designed around short (~150bp) paired-end reads, where mate-pair information helps disambiguate repeats and adapter read-through rarely reaches deep enough to matter. **SE600 data is a different regime**: single-end, ~600bp reads with no mate to cross-check against, long enough to routinely read past a shorter physical insert into ligation adapter, sample index, and even the barcode itself, and assembled at a scale (millions of barcodes per run) where noise that would be a rare edge case at low throughput becomes a routine, high-volume failure mode. Two problems fall directly out of this:

- **De Bruijn graph assembly requires k-mer coverage ≥ 2.** Low-depth barcodes (a common, non-negligible fraction of any real barcode-partitioned dataset) simply cannot be assembled — the graph fragments or never forms.
- **Forking a full assembler subprocess per barcode does not scale.** At 1.5–3M barcodes per run, subprocess-spawn overhead alone becomes the dominant cost, independent of the actual assembly work.

`denovo_OLC` addresses both by implementing the classic OLC paradigm (overlap detection → layout → consensus) directly in-process, in pure Python, with no external assembler dependency and no per-barcode subprocess fork — including graceful assembly from as little as a **single read** — while treating SE600-specific noise (boundary noise, combined independent read errors, adapter/index/barcode read-through) as first-class problems to detect and handle rather than edge cases to ignore.

## Core innovations

### 1. Output-sensitive candidate discovery via boundary k-mer inverted indices
Rather than re-scanning the full read pool for every contig-extension step (the naive OLC approach, which degrades to O(attempts × pool size) on barcodes with poor pairwise mergeability), a per-UMI inverted index maps boundary k-mers to candidate reads once, and every extension step queries it directly. On a documented worst-case real-world barcode, this cut single-UMI assembly time from **100.2s to 5.07s (11.7×)** while producing a byte-for-byte identical assembly (verified by digest equality against the unindexed baseline across isolated and 500-barcode batch tests) — the index only prunes the search, it never changes which overlaps are accepted.

### 2. Bidirectional internal k-mer anchoring
Real reads carry noise at their literal 5′/3′ ends (adapter read-through, sequencing error enrichment near read termini). A naive OLC that only tries to match read *boundaries* silently fails to extend past such noise even when a clean, unambiguous overlap exists a few bases further in. `denovo_OLC` adds an internal-anchor fallback — indexed in both the forward and reverse direction — that locates a verified internal anchor and extends from there, recovering merges that boundary-only matching would miss entirely.

### 3. Collective pileup rescue with cross-attempt evidence sharing
Pairwise, single-read overlap comparison has a structural weakness: two reads that each carry independent sequencing errors combine both error rates when compared directly, often pushing a genuinely valid overlap past the mismatch-rate threshold. Rather than relaxing that threshold globally (which would open the door to chimeric misassembly), `denovo_OLC` falls back — only when ordinary pairwise extension stalls — to anchoring multiple reads simultaneously via long, low-collision k-mers and extending only the region with pileup (multi-read) support. Evidence can be shared across separate greedy build attempts on the same barcode, not just within one, letting reads consumed by an earlier attempt still contribute supporting votes to a later one. A progress invariant (every accepted extension must consume at least one previously-unused read) guarantees termination even under adversarial repetitive input.

### 4. Lightweight majority-vote consensus polish
After assembly, a k-mer-offset voting pass re-aligns every input read against the draft contig (no full dynamic-programming realignment) and flips a position only when there is majority, quorum-backed disagreement — correcting substitution errors introduced during greedy extension without the cost or complexity of a full multiple-sequence-alignment consensus step.

### 5. Deterministic-by-construction output
Python's per-process string-hash randomization means naive use of `set()` for deduplication can silently make seed selection — and therefore the entire greedy assembly trajectory — non-reproducible across runs of *identical* input. Every ordering-sensitive step in `denovo_OLC` uses an explicit, content-based deterministic tie-break, so the same input always produces the same output regardless of interpreter hash seed, process start method, or worker count.

### 6. Adapter/index/barcode read-through detection
This is where SE600's read length becomes a double-edged sword: a 600bp single-end read routinely runs past a shorter physical insert, straight through ligation adapter and sample index, and — in rolling-circle-amplified (DNB-style) library preps — back into the barcode itself. Left untrimmed, this technical read-through is assembled as if it were biological sequence, inflating apparent contig length and fragmenting what should be a single contiguous assembly. `denovo_OLC`'s validation tooling includes reference-free forensic methods (positional-invariance testing, abundance-matched conserved-region controls, and cross-referencing against reference sequence databases) to distinguish genuine adapter contamination from true biological conserved regions — critical for correctly attributing this SE600-specific class of length/fragmentation artifact rather than misdiagnosing it as an assembly bug.

### 7. Production-grade multiprocessing correctness
Handles the divergent semantics of `fork` (Linux default) vs. `spawn` (macOS/Windows default) multiprocessing start methods correctly — configuration state is explicitly propagated to every worker rather than relying on `fork`'s copy-on-write inheritance, which silently breaks under `spawn`. Task scheduling is tuned for the heterogeneous per-barcode cost distribution inherent to biological data (most barcodes cost milliseconds; a small tail costs orders of magnitude more), avoiding the worker starvation that a naive fixed-chunksize scheduler produces under that distribution.

## Design philosophy

- **No external assembler dependency.** Pure standard library — no compiled extensions, no subprocess management, no non-deterministic third-party graph library behavior to reason about.
- **Correctness has priority over aggressive heuristics.** Every fallback mechanism (internal anchoring, collective rescue) is deliberately ordered *after* the safest, most conservative option succeeds or is exhausted, and every optimization is validated against a fixed-output-digest regression before being accepted.
- **Validated on real data, not just synthetic benchmarks.** The synthetic test suite (unit tests covering overlap detection, chimeric-merge safety, deterministic tie-breaking, multiprocessing config propagation, and more) is paired with forensic analysis against real sequencing data to catch failure modes synthetic tests structurally cannot — asymmetric read-length distributions, platform-specific adapter constructs, and genuinely pathological real-world depth outliers.

## Usage

```bash
python3 src/denovo_seed_olc.py \
  --sequence_type se \
  --num_processes 30 \
  --min_ctg_len 400 \
  --r2 data_R2_sgrep.tsv
```

See the module docstrings for the full configuration surface (overlap/mismatch thresholds, internal-anchor and collective-rescue toggles, polish parameters, output caps).

## Status

Actively used in production for per-UMI 16S rRNA amplicon and metagenomic fragment reconstruction at multi-million-UMI scale.

## License

Released under [CC BY-NC 4.0](LICENSE) (Attribution-NonCommercial) — free to share and adapt with attribution; commercial use is not permitted. See [LICENSE](LICENSE) for the full text.
