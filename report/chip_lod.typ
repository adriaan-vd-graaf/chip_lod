// chip_lod report. Every number, table and figure is read from results/ — recompile after a rerun.
// Compile from the repo root: typst compile --root . report/chip_lod.typ

#let res = "../results/python/"
#let s = json(res + "summary.json")
#let hc = json("../results/hash_check.json")
#let F = s.f

// ---------- helpers ----------

#let load(name) = {
  let rows = csv(res + name)
  let h = rows.first()
  rows.slice(1).map(r => h.zip(r).to-dict())
}
#let near(a, b) = calc.abs(float(a) - float(b)) <= 1e-9 * calc.max(1.0, calc.abs(float(b)))
#let isna(x) = x == none or x == "NA"
#let fmt(x, d: 2) = if isna(x) { [--] } else { str(calc.round(float(x), digits: d)) }
#let pct(x, d: 2) = if isna(x) { [--] } else { str(calc.round(float(x) * 100, digits: d)) + "%" }
#let istr(x) = if isna(x) { [--] } else { str(calc.round(float(x))) }
#let sci(x, d: 1) = {
  let v = float(x)
  if v == 0 { return [0] }
  let e = calc.floor(calc.log(v, base: 10))
  let m = v / calc.pow(10.0, e)
  [#str(calc.round(m, digits: d)) × 10#super[#str(e)]]
}
#let prob(x) = {
  if isna(x) { return [--] }
  let v = float(x)
  if v == 0 { [0] } else if v == 1 { [1] } else if v < 0.001 { sci(v) } else if v > 0.9999 { [> 0.9999] } else { str(calc.round(v, digits: 4)) }
}
#let ppl(x) = { // a "number of people" in the natural-frequency examples
  let v = float(x)
  if v >= 10 { str(calc.round(v, digits: 0)) } else if v >= 0.01 { str(calc.round(v, digits: 1)) } else { sci(v) }
}
#let dfmt(x) = if x == none { [never (within the grid)] } else { [#x reads] }
#let keybox(body) = block(fill: luma(240), inset: 10pt, radius: 4pt, width: 100%, body)
#let note(body) = block(inset: (left: 8pt), stroke: (left: 2pt + luma(180)), body)

#set document(title: "chip_lod: how many reads make a CHIP call?", author: "chip_lod")
#set page(paper: "a4", margin: 2.2cm, numbering: "1")
#set text(size: 10.5pt)
#set par(justify: true)
#set heading(numbering: "1.")
#show table: set text(size: 9pt)
#show figure.caption: set text(size: 9pt)

#align(center)[
  #text(size: 18pt, weight: "bold")[chip_lod: how many alternate reads make a CHIP call?] \
  #v(2pt)
  #text(size: 11pt)[A plain-language guide to the model, the analyses and how to read them]
]

#v(6pt)
#keybox[
  *The question.* At one known position (for example a _DNMT3A_ or _JAK2_ hotspot), a sequencer reports $n$ reads,
  and $k$ of them show the variant base. Some of those reads can be sequencing errors. How many variant reads do we
  need before we can trust the call, and how small a variant can we detect at all?

  *The answer in short.* With #F.reference_depth reads and an error rate of #F.reference_error_rate, a call needs at least
  *#F.k_star variant reads*. With that rule a variant at *#pct(F.lod_vaf) VAF* or more is found in #pct(F.power, d: 0) of
  samples: that is the *limit of detection* (LoD). The LoD rises with the error rate (from #pct(F.by_error_rate.first().lod_vaf)
  at $e$ = #F.by_error_rate.first().error_rate to #pct(F.by_error_rate.last().lod_vaf) at $e$ = #F.by_error_rate.last().error_rate)
  and falls with depth. Detecting a #pct(F.required_depth.last().target_vaf, d: 0) variant needs
  #F.required_depth.last().min_depth reads, and a #pct(F.required_depth.first().target_vaf, d: 1) variant
  #F.required_depth.first().min_depth reads.
]

#outline(depth: 1, indent: 1em)

= The question

Clonal haematopoiesis of indeterminate potential (CHIP) is defined by a somatic mutation in blood at a variant
allele fraction (VAF) of at least #pct(s.model.chip_vaf, d: 0). Blood cells are diploid and the variant sits on one
of the two copies, so a VAF of 2% means about 4% of cells carry it.

Sequencing is not perfect. Even at a position where nobody carries the variant, a small fraction of reads show
the variant base because of sequencing or PCR errors. With an error rate $e$ = #F.reference_error_rate per base,
and a third of errors landing on any particular wrong base, a site with #F.reference_depth reads is expected to show
*#fmt(F.lambda_bg, d: 3) error reads* on average. So when we see 1, 2 or 5 variant reads, is that a variant or noise?

This report answers it with classical (frequentist) statistics only: a significance test against the error
background, and the limit of detection that follows from it. We deliberately do not use a Bayesian analysis,
because that needs prior information (how common the variant is, and at what VAF) that we do not have.

All results come from the `chip_lod` package (Python, with an R port that produces byte-identical files). All
settings are in `config/params.json`.

#let f1 = load("f1_lod_grid.csv")
#let f2 = load("f2_required_depth.csv")
#let f3 = load("f3_power_curve.csv")
#let f4 = load("f4_sim_power.csv")
#let f5 = load("f5_lod_function_example.csv")
#let f-errs = f1.map(r => r.error_rate).dedup()
#let f-depths = f1.map(r => r.depth).dedup()
#let ref-ex = f5.filter(r => int(r.n) == F.reference_depth)

= How the frequentist test works

The frequentist approach asks one question: *if there were no variant, how surprising would these reads be?*
It needs only a model of the errors. Here is the procedure at depth #F.reference_depth with $e$ = #F.reference_error_rate.

+ *Background.* Errors produce on average $lambda_0 = n dot e slash 3$ = #fmt(F.lambda_bg, d: 3) variant reads. The
  actual number varies from sample to sample and follows a Poisson distribution, the standard model for counts
  of rare, independent events.
+ *p-value.* The p-value of $k$ reads is the probability that errors alone produce $k$ _or more_ variant reads.
  #for r in ref-ex [With $k$ = #r.k the p-value is #prob(r.pvalue). ]
+ *Critical value $k^*$.* We call the variant when the p-value is at most $alpha$ = #F.alpha, i.e. when errors would
  produce that many reads in fewer than #calc.round(F.alpha * 100)% of variant-free samples. The smallest such
  count is $k^*$ = *#F.k_star*. Because read counts are whole numbers, the actual false-positive rate is
  #prob(F.actual_alpha), a little below $alpha$.
+ *Power.* A real variant at VAF $p$ produces variant reads at the rate $r(p) = p(1-e) + (1-p) e slash 3$ per read
  (correct reads from the mutated copy plus misread normal reads). The *power* is the chance that such a variant
  reaches $k^*$ reads and is called.
+ *Limit of detection.* The LoD is the smallest VAF that is called in at least #pct(F.power, d: 0) of samples.
  At depth #F.reference_depth this is *#pct(F.lod_vaf)*.

#note[
  *How to read this.* $k^*$ is a property of the depth and the error rate, and is fixed before looking at the data.
  The LoD is a property of the assay: it says which variants the assay reliably detects, not whether this
  particular sample carries one. A variant below the LoD can still be called, just not reliably; a sample
  without a call can still carry a variant below the LoD.
]

The next sections vary the three ingredients one at a time (error rate, depth and VAF), focusing on VAFs up to
#pct(F.vaf_max, d: 0), the range where calling is difficult. Throughout, $alpha$ = #F.alpha and the required power is
#pct(F.power, d: 0).

= Influence of the error rate

#figure(
  table(
    columns: 5,
    align: right,
    table.header([error rate $e$], [background $lambda_0$], [$k^*$], [actual false-positive rate], [LoD VAF]),
    ..F.by_error_rate.map(r => ([#r.error_rate], fmt(r.lambda_bg, d: 3), [#r.k_star], prob(r.actual_alpha),
      pct(r.lod_vaf))).flatten(),
  ),
  caption: [The effect of the error rate at depth #F.reference_depth (`f1_lod_grid.csv`).],
)

A higher error rate means more background reads, so a call needs more variant reads ($k^*$ rises from
#F.by_error_rate.first().k_star to #F.by_error_rate.last().k_star), and the LoD rises with it. Going from the
cleanest assay ($e$ = #F.by_error_rate.first().error_rate, typical of UMI/duplex consensus reads) to the noisiest
($e$ = #F.by_error_rate.last().error_rate) raises the LoD about
#calc.round(float(F.by_error_rate.last().lod_vaf) / float(F.by_error_rate.first().lod_vaf), digits: 1)-fold.

The "actual false-positive rate" column shows the effect of whole-number thresholds. The test promises at most
#pct(F.alpha, d: 0) false positives, but because $k^*$ jumps in whole steps, the real rate is anywhere between almost
0 and #pct(F.alpha, d: 0).

#figure(
  table(
    columns: 1 + f-errs.len(),
    align: right,
    table.header([depth], ..f-errs.map(e => [$e$ = #fmt(e, d: 4)])),
    ..f-depths.map(d => ([#d], ..f-errs.map(e => {
      let r = f1.find(r => r.depth == d and r.error_rate == e)
      [#pct(r.lod_vaf) (#r.k_star)]
    }))).flatten(),
  ),
  caption: [LoD VAF for every combination of depth and error rate, with $k^*$ in brackets (`f1_lod_grid.csv`).],
)

= Influence of the read depth

#figure(
  image("../results/figures/fig_f1_lod_depth.svg", width: 72%),
  caption: [Limit of detection against depth for six error rates (log–log). The dotted line is the 2% CHIP
    threshold.],
)

More reads make the variant reads stand out from the background, so the LoD falls with depth. While errors are
negligible ($k^*$ stays at 1 or 2), it falls almost in proportion to $1 slash n$: ten times the depth gives a ten
times smaller LoD. Once the background amounts to several reads, the variant has to rise above the random scatter
of those error reads, and the LoD falls more slowly, closer to $1 slash sqrt(n)$. At $e$ = #F.reference_error_rate the LoD goes from
#pct(f1.find(r => near(r.error_rate, F.reference_error_rate)).lod_vaf) at #f-depths.first() reads to
#pct(f1.filter(r => near(r.error_rate, F.reference_error_rate)).last().lod_vaf) at #f-depths.last() reads.
The lines are not perfectly smooth: each time the depth grows enough for $k^*$ to step up by one read, the LoD
briefly gets _worse_ before improving again. This is the "saw-tooth" effect of whole-number thresholds.

*How deep do I need to sequence?* The table below inverts the question. It gives the smallest depth at which a
variant of a given VAF is detected with #pct(F.power, d: 0) power.

#let f2-targets = f2.map(r => r.target_vaf).dedup()
#figure(
  table(
    columns: 1 + f2-targets.len(),
    align: right,
    table.header([error rate $e$], ..f2-targets.map(v => [VAF #pct(v, d: 1)])),
    ..f-errs.map(e => ([#fmt(e, d: 4)], ..f2-targets.map(v => {
      let r = f2.find(r => r.error_rate == e and r.target_vaf == v)
      [#istr(r.min_depth) ($k^*$ = #istr(r.k_star))]
    }))).flatten(),
  ),
  caption: [Smallest depth with #pct(F.power, d: 0) power to detect each VAF (`f2_required_depth.csv`).],
)

#figure(
  image("../results/figures/fig_f3_required_depth.svg", width: 66%),
  caption: [Depth needed for #pct(F.power, d: 0) power against error rate, for three target VAFs.],
)

#let nonmono = f2-targets.filter(v => {
  let ds = f-errs.map(e => f2.find(r => r.error_rate == e and r.target_vaf == v).min_depth).map(float)
  range(ds.len() - 1).any(i => ds.at(i + 1) < ds.at(i))
})
Small VAFs are expensive: halving the VAF roughly doubles (or more) the depth needed, and a noisy assay can
need several times more reads than a clean one.
#if nonmono.len() > 0 [
  One result looks paradoxical: for #nonmono.map(v => pct(v, d: 1)).join(", ", last: " and ") VAF, a
  slightly _higher_ error rate sometimes needs _fewer_ reads. This happens when both error rates use the same
  $k^*$. The extra error reads then also count towards the variant's total and help it cross the threshold.
  The price is a higher actual false-positive rate (previous section), so the noisier assay is not really better.
]

= Influence of the VAF

#figure(
  image("../results/figures/fig_f2_power.svg", width: 90%),
  caption: [Power (probability of a call) against true VAF up to #pct(F.vaf_max, d: 0). Left: depths at $e$ =
    #F.reference_error_rate. Right: error rates at depth #F.reference_depth. The dotted line is #pct(F.power, d: 0)
    power (`f3_power_curve.csv`).],
)

#let f3-ref = f3.filter(r => near(r.error_rate, F.reference_error_rate))
#let f3-depths = f3-ref.map(r => r.depth).dedup()
#let show-vafs = (0.001, 0.0025, 0.005, 0.01, 0.02).filter(v => v <= float(F.vaf_max) + 1e-12)
#figure(
  table(
    columns: 1 + f3-depths.len(),
    align: right,
    table.header([true VAF], ..f3-depths.map(d => [n = #d ($k^*$ = #f3-ref.find(r => r.depth == d).k_star)])),
    ..show-vafs.map(v => ([#pct(v, d: 2)], ..f3-depths.map(d => {
      let r = f3-ref.find(r => r.depth == d and near(r.true_vaf, v))
      if r == none [--] else [#pct(r.power, d: 1)]
    }))).flatten(),
  ),
  caption: [Probability of calling a variant of a given VAF at $e$ = #F.reference_error_rate (selected rows of
    `f3_power_curve.csv`).],
)

#let f1-ref = f1.filter(r => near(r.error_rate, F.reference_error_rate))
#let d01 = f1-ref.find(r => float(r.lod_vaf) <= 0.001)
For a fixed assay, power rises steeply with VAF: below the LoD a variant is called only some of the time, and above it
almost always. Between about 0.1% and 1% VAF, the region of interest for early clonal haematopoiesis, the answer
depends strongly on depth. At $e$ = #F.reference_error_rate, a variant at the
#pct(F.required_depth.last().target_vaf, d: 0) CHIP threshold is reliably detected from
#F.required_depth.last().min_depth reads onwards, #pct(F.required_depth.first().target_vaf, d: 1) needs
#F.required_depth.first().min_depth reads, and 0.1% #if d01 == none [is not reached at any depth in our grid
(up to #f-depths.last() reads)] else [is first reached at #d01.depth reads in our depth grid].

= A function for the frequentist LoD

`lod_frequentist(n, k, e, alpha, power)` applies the whole analysis to a list of sites. Given vectors of depths
$n$ and variant-read counts $k$, it returns for every site: the background $lambda_0$, the p-value, $k^*$, whether the
site is called, and the LoD at that depth.

#grid(columns: 2, gutter: 10pt,
  ```python
  from chip_lod.lod import lod_frequentist
  res = lod_frequentist(n=[500, 1000, 2000],
                        k=[3, 1, 3], e=0.001)
  res["called"], res["lod_vaf"]
  ```,
  ```r
  source("r/load.R")
  lod_frequentist(n = c(500, 1000, 2000),
                  k = c(3, 1, 3), e = 0.001)
  # data.frame: n, k, lambda_bg, pvalue,
  #             k_star, called, lod_vaf
  ```,
)

#figure(
  table(
    columns: 7,
    align: right,
    table.header([$n$], [$k$], [$lambda_0$], [p-value], [$k^*$], [called], [LoD VAF]),
    ..f5.map(r => ([#r.n], [#r.k], fmt(r.lambda_bg, d: 3), prob(r.pvalue), [#r.k_star],
      if r.called == "1" [yes] else [no], pct(r.lod_vaf))).flatten(),
  ),
  caption: [Output of `lod_frequentist` for example sites at $e$ = #F.reference_error_rate (`f5_lod_function_example.csv`).],
)

*Methods paragraph* (for a paper; also in `spec/methods_lod_frequentist.txt`):

#block(inset: 10pt, stroke: 0.5pt + luma(160), radius: 4pt, width: 100%, text(size: 9.5pt,
  read("../spec/methods_lod_frequentist.txt")))

= In practice: a panel at depth #s.a3.depth

#let a3r = load("a3_read_table.csv")
#let a3t = load("a3_thresholds.csv")
#let a3s = load("a3_sim_sensitivity.csv")
#let sim3-depths = csv(res + "sim_a3.csv").slice(1).map(r => int(r.at(2)))
#let deep = a3t.filter(r => near(r.alpha, s.a3.alpha)).last()
#let a3-alphas = a3t.map(r => r.alpha).dedup()
#let a3-depths = a3t.map(r => r.depth).dedup()
#let ks-strict = a3t.find(r => int(r.depth) == s.a3.depth and not near(r.alpha, s.a3.alpha))

A typical targeted panel gives about #s.a3.depth reads at a hotspot, and the depth varies from sample to sample.
At $e$ = #s.a3.error_rate the background is #fmt(s.a3.lambda_bg, d: 3) error reads.

#figure(
  table(
    columns: 4,
    align: right,
    table.header([variant reads $k$], [p-value], [called at $alpha$ = #fmt(a3-alphas.first(), d: 3)?],
      [called at $alpha$ = #fmt(a3-alphas.last(), d: 3)?]),
    ..a3r.filter(r => float(r.alt_reads) <= 8).map(r => ([#r.alt_reads], prob(r.pvalue),
      if float(r.pvalue) <= float(a3-alphas.first()) [yes] else [no],
      if float(r.pvalue) <= float(a3-alphas.last()) [yes] else [no])).flatten(),
  ),
  caption: [p-values for 0–8 variant reads at depth #s.a3.depth (`a3_read_table.csv`, which runs to $k$ =
    #a3r.last().alt_reads).],
)

#figure(
  table(
    columns: 1 + a3-alphas.len() * 3,
    align: right,
    table.header(
      table.cell(rowspan: 2)[depth],
      ..a3-alphas.map(a => table.cell(colspan: 3, align: center)[$alpha$ = #fmt(a, d: 3)]),
      ..a3-alphas.map(a => ([$k^*$], [actual FP rate], [LoD VAF])).flatten()),
    ..a3-depths.map(d => ([#d], ..a3-alphas.map(a => {
      let r = a3t.find(r => r.depth == d and r.alpha == a)
      ([#r.k_star], prob(r.actual_alpha), pct(r.lod_vaf))
    }).flatten())).flatten(),
  ),
  caption: [Critical value, actual false-positive rate and LoD around depth #s.a3.depth for a lenient and a strict
    $alpha$ (`a3_thresholds.csv`).],
)

A stricter significance level protects better against false calls but costs sensitivity: at depth #s.a3.depth,
$alpha$ = #fmt(ks-strict.alpha, d: 3) needs #ks-strict.k_star reads instead of #s.a3.k_star and raises the LoD from
#pct(s.a3.lod_vaf) to #pct(ks-strict.lod_vaf). The LoD does not always fall as depth rises. When $k^*$ steps up by
one read, the LoD can jump up slightly#if float(deep.lod_vaf) > float(s.a3.lod_vaf) [ (here from #pct(s.a3.lod_vaf) at
depth #s.a3.depth to #pct(deep.lod_vaf) at depth #deep.depth, $alpha$ = #s.a3.alpha)].

*Simulated panel.* We simulated sites whose depth varies uniformly between #calc.min(..sim3-depths) and
#calc.max(..sim3-depths), and called each one with its own $k^*$ ($alpha$ = #s.a3.alpha).

#figure(
  table(
    columns: 5,
    align: right,
    table.header([true VAF], [sites], [mean depth], [call rate], [analytic power at #s.a3.depth]),
    ..a3s.map(r => (pct(r.true_vaf), [#r.n_sites], fmt(r.mean_depth, d: 0), pct(r.call_rate_freq, d: 1),
      pct(r.analytic_power_at_1000, d: 1))).flatten(),
  ),
  caption: [Call rates in simulated data with variable depth (`a3_sim_sensitivity.csv`). The VAF 0 row is the
    false-positive rate.],
)

#figure(
  image("../results/figures/fig_a3_sensitivity.svg", width: 70%),
  caption: [Call rate against true VAF for sites with variable depth (99% Wilson intervals), with the analytic power
    at depth #s.a3.depth.],
)

With no variant, #pct(s.a3.sim_fp_rate_freq, d: 1) of sites were called, #if float(s.a3.sim_fp_rate_freq) <= float(s.a3.alpha) [below] else [*above*] $alpha$ = #s.a3.alpha. The
analytic power assumes exactly #s.a3.depth reads. Deeper simulated sites need more reads ($k^*$ = #deep.k_star at
depth #deep.depth), which is why the simulated call rates sit slightly below it: when depth varies between
samples, the LoD of the panel is set by its deeper sites as well as its average.

= Reproducibility

The package has two independent implementations, Python and R, written to produce *byte-identical* result files.
Both use the same hand-written random number generator, explicit left-to-right sums and the same number
formatting. `scripts/compare_outputs.sh` compares the SHA-256 hash of every file.

#let hfiles = hc.files.keys()
#let n-match = hfiles.filter(f => hc.files.at(f).match).len()
#block(fill: if hc.all_match { rgb("#e6f2e6") } else { rgb("#f8e0e0") }, inset: 8pt, radius: 4pt, width: 100%)[
  *#if hc.all_match [Python and R agree:] else [Python and R DISAGREE:]* #n-match of #hfiles.len() result files
  have identical SHA-256 hashes.
]

#table(
  columns: 3,
  table.header([file], [SHA-256 (first 16 hex digits)], [match]),
  ..hfiles.map(f => {
    let e = hc.files.at(f)
    (raw(f), raw(if e.python == none { "missing" } else { e.python.slice(0, 16) }), if e.match [yes] else [*no*])
  }).flatten(),
)

Every number in this report is read from `results/python/` when the document is compiled, so rerunning `make all`
updates everything. The model is specified in `spec/MODEL.md`, and every judgement call is listed in
`spec/DECISIONS.md`.

= Limitations

- *Poisson approximation.* Read counts are really binomial, and the Poisson is slightly too wide. The simulations
  show the difference is negligible at these depths and VAFs.
- *No overdispersion.* Real error counts vary more from sample to sample than a Poisson allows (batch, library and
  position effects). A beta-binomial or negative-binomial background would give larger, more honest $k^*$ values.
- *The error rate is site- and context-specific.* A single $e$ = #F.reference_error_rate is a placeholder. Error rates
  depend on the base change, the sequence context and the assay; the error-rate section shows how strongly this matters. In practice,
  estimate the background at each hotspot from a *panel of normals*.
- *PCR duplicates and UMIs.* The model treats reads as independent. Duplicates from the same molecule are not
  independent and make a variant look better supported than it is. Count unique molecules (UMI families), not raw
  reads; UMI consensus also lowers $e$ considerably.
- *Diploid, heterozygous assumption.* VAF = cell fraction / 2 assumes a heterozygous variant on a normal copy number.
  Loss of heterozygosity, copy-number changes or a homozygous variant break this link.
